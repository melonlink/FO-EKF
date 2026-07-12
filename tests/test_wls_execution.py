import copy
from fractions import Fraction

import numpy as np
import pytest

from fo_ekf.wls_execution import (
    bounded_sample_indices,
    build_wls_execution_record,
    canonical_hash,
    replay_wls_execution_record,
    validate_execution_manifest_binding,
)


class _FloatSubclass(float):
    pass


def _execution() -> tuple[dict, np.ndarray]:
    beats = tuple(range(100, 261, 10))
    samples = np.arange(beats[0], beats[-1])
    phase = 2.0 * np.pi * (samples - beats[0]) / 10.0
    physical = 2.0 + 0.7 * np.cos(phase) - 0.2 * np.sin(phase)
    gain = 1000.0
    digital = np.rint(physical * gain).astype(np.int64)
    record = build_wls_execution_record(
        digital,
        beats,
        sampling_frequency=250.0,
        adc_gain=gain,
        baseline=0,
        harmonics=(1, 2, 3),
        protocol_sha256="1" * 64,
        record_sha256="2" * 64,
        window_sha256="3" * 64,
    )
    return record, digital


def test_execution_record_binds_payload_statistics_center_and_ball_bound() -> None:
    record, digital = _execution()

    assert record["binding"]["protocol_sha256"] == "1" * 64
    assert record["weights"]["raw_definition"] == "identity raw weights w_raw_n=1"
    assert record["window"]["retained_harmonics"] == [-3, -2, -1, 0, 1, 2, 3]
    assert len(record["sufficient_statistics"]["G"]) == 7
    assert len(record["sufficient_statistics"]["b"]) == 7
    assert len(record["result"]["C_tilde"]) == 7
    assert len(record["result"]["z_tilde"]) == 3
    assert record["result"]["z_tilde"][0]["harmonic"] == 1
    assert np.isclose(record["result"]["z_tilde"][0]["real"], 0.35, atol=1.0e-3)
    assert np.isclose(record["result"]["z_tilde"][0]["imag"], 0.1, atol=1.0e-3)
    assert record["result"]["joint_solve_rounding_enclosure"]["status"] == "ENCLOSED"
    rounding = record["result"]["joint_solve_rounding_enclosure"]
    exact_upper = rounding["maximum_absolute_error_upper_bound_dyadic"]
    assert Fraction(
        int(exact_upper["numerator"]), int(exact_upper["denominator"])
    ) == Fraction.from_float(rounding["maximum_absolute_error_upper_bound"])
    assert set(record["environment"]) == {
        "python",
        "numpy",
        "python_flint",
        "flint",
    }
    assert replay_wls_execution_record(record, digital)["status"] == "MATCH"


def test_replay_detects_sample_tamper_and_document_tamper() -> None:
    record, digital = _execution()
    changed = digital.copy()
    changed[5] += 1

    replay = replay_wls_execution_record(record, changed)
    assert replay["status"] == "MISMATCH"
    assert "payload_hashes" in replay["mismatched_fields"]

    tampered = copy.deepcopy(record)
    tampered["sufficient_statistics"]["b"][0]["real"] += 0.01
    replay = replay_wls_execution_record(tampered, digital)
    assert replay == {"status": "MISMATCH", "reason": "stored_execution_hash_mismatch"}


def test_replay_rejects_invalid_hash_binding_fail_closed() -> None:
    record, digital = _execution()
    tampered = copy.deepcopy(record)
    tampered["binding"]["window_sha256"] = "not-a-hash"
    tampered["execution_sha256"] = "0" * 64

    replay = replay_wls_execution_record(tampered, digital)

    assert replay["status"] == "MISMATCH"


def test_resource_and_exact_type_guards_fail_closed() -> None:
    record, digital = _execution()
    arguments = {
        "sampling_frequency": 250.0,
        "adc_gain": 1000.0,
        "baseline": 0,
        "harmonics": (1, 2, 3),
        "protocol_sha256": "1" * 64,
        "record_sha256": "2" * 64,
        "window_sha256": "3" * 64,
    }
    with pytest.raises(ValueError, match="nonempty"):
        build_wls_execution_record(digital[:0], (), **arguments)
    with pytest.raises(ValueError, match="exact integers"):
        build_wls_execution_record(
            digital, tuple(range(100, 261, 10)), **(arguments | {"harmonics": (True,)})
        )
    with pytest.raises(ValueError, match="exact integers"):
        build_wls_execution_record(
            digital,
            tuple(range(100, 261, 10)),
            **(arguments | {"harmonics": (1.0,)}),
        )
    with pytest.raises(ValueError, match="ordered and unique"):
        build_wls_execution_record(
            digital, tuple(range(100, 261, 10)), **(arguments | {"harmonics": (2, 1)})
        )
    with pytest.raises(ValueError, match="baseline"):
        build_wls_execution_record(
            digital, tuple(range(100, 261, 10)), **(arguments | {"baseline": True})
        )
    with pytest.raises(ValueError, match="integer dtype"):
        build_wls_execution_record(
            digital.astype(np.float64), tuple(range(100, 261, 10)), **arguments
        )
    with pytest.raises(ValueError, match="numeric and not boolean"):
        build_wls_execution_record(
            digital,
            tuple(range(100, 261, 10)),
            **arguments,
            front_end_response_by_harmonic={1: True, 2: 1.0, 3: 1.0},
        )

    long_beats = tuple(range(0, 5001, 312))[:16] + (5000,)
    selected = bounded_sample_indices(0, 5000)
    long_record = build_wls_execution_record(
        np.zeros(selected.size, dtype=np.int64), long_beats, **arguments
    )
    assert long_record["window"]["sample_count"] == 4096
    assert long_record["window"]["full_span_sample_count"] == 5000
    assert (
        replay_wls_execution_record(long_record, np.zeros(selected.size, dtype=np.int64))["status"]
        == "MATCH"
    )
    assert record["window"]["sample_count"] < 4096


def test_replay_handles_resigned_huge_beat_and_deep_tree_without_crash() -> None:
    record, digital = _execution()
    huge = copy.deepcopy(record)
    huge["window"]["beat_samples"][-1] = 2**100
    unsigned = dict(huge)
    unsigned.pop("execution_sha256")
    huge["execution_sha256"] = canonical_hash(unsigned)

    replay = replay_wls_execution_record(huge, digital)
    assert replay["status"] == "UNKNOWN"
    assert "ValueError" in replay["reason"]

    deep = copy.deepcopy(record)
    cursor = deep
    for _ in range(40):
        cursor["nested"] = {}
        cursor = cursor["nested"]
    replay = replay_wls_execution_record(deep, digital)
    assert replay == {"status": "UNKNOWN", "reason": "document_depth_budget_exceeded"}

    numeric_subclass = copy.deepcopy(record)
    numeric_subclass["window"]["adc_gain"] = _FloatSubclass(numeric_subclass["window"]["adc_gain"])
    replay = replay_wls_execution_record(numeric_subclass, digital)
    assert replay == {"status": "UNKNOWN", "reason": "document_non_plain_json_value"}


def test_manifest_source_execution_cross_binding_is_one_record_exact() -> None:
    execution, _digital = _execution()
    wrapper = {
        "record": "synthetic",
        "role": "test",
        "split": "validation",
        "window_index": 7,
        "execution": execution,
    }
    window = execution["window"]
    manifest = {
        "record": "synthetic",
        "role": "test",
        "split": "validation",
        "window_index": "7",
        "window_sha256": execution["binding"]["window_sha256"],
        "wls_execution_sha256": execution["execution_sha256"],
        "start_sample": str(window["start_sample_inclusive"]),
        "end_sample_exclusive": str(window["end_sample_exclusive"]),
        "sample_count": str(window["sample_count"]),
        "full_span_sample_count": str(window["full_span_sample_count"]),
    }

    assert (
        validate_execution_manifest_binding(
            wrapper,
            manifest,
            external_record_sha256="2" * 64,
            protocol_sha256="1" * 64,
        )
        == []
    )
    manifest["wls_execution_sha256"] = "0" * 64
    failures = validate_execution_manifest_binding(
        wrapper,
        manifest,
        external_record_sha256="f" * 64,
        protocol_sha256="1" * 64,
    )
    assert failures == ["binding_record_sha256", "execution_sha256"]
