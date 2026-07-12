import copy
from fractions import Fraction

import numpy as np
import pytest
from flint import ctx

from fo_ekf.wls_execution import build_wls_execution_record, canonical_hash
from fo_ekf.wls_interval_certificate import (
    _validated_int64_payload,
    build_wls_interval_certificate,
    certified_radius_floors,
    replay_wls_interval_certificate,
)


class _RatioZero(float):
    def as_integer_ratio(self) -> tuple[int, int]:
        return (0, 1)


def _execution(
    *,
    front_end: dict[int, complex] | None = None,
) -> tuple[dict, np.ndarray]:
    beats = tuple(range(100, 261, 10))
    samples = np.arange(beats[0], beats[-1])
    phase = 2.0 * np.pi * (samples - beats[0]) / 10.0
    physical = 2.0 + 0.7 * np.cos(phase) - 0.2 * np.sin(phase)
    gain = 1000.0
    digital = np.rint(physical * gain).astype(np.int64)
    execution = build_wls_execution_record(
        digital,
        beats,
        sampling_frequency=250.0,
        adc_gain=gain,
        baseline=0,
        harmonics=(1, 2, 3),
        protocol_sha256="1" * 64,
        record_sha256="2" * 64,
        window_sha256="3" * 64,
        front_end_response_by_harmonic=front_end,
    )
    return execution, digital


def test_exact_payload_certificate_closes_full_numeric_wls_path() -> None:
    execution, digital = _execution()

    certificate = build_wls_interval_certificate(execution, digital, precision_bits=192)

    assert certificate["relation"] == "CERTIFIED_INTERVAL"
    assert certificate["binding"]["execution_sha256"] == execution["execution_sha256"]
    assert (
        certificate["binding"]["digital_sample_sha256"]
        == execution["payload_hashes"]["digital_sample_sha256"]
    )
    assert certificate["numeric_proof"]["determinant_zero_excluded"] is True
    assert certificate["mathematical_contract"]["retained_harmonics"] == [
        -3,
        -2,
        -1,
        0,
        1,
        2,
        3,
    ]
    assert "G_and_b_direct_accumulation_rounding" in certificate["coverage"]
    assert "G_and_b_accumulation_rounding" not in certificate["exclusions"]
    assert len(certificate["response_disks"]) == 3
    assert (
        "not asserted to be the minimal true error"
        in certificate["response_disks"][0]["radius_floor_semantics"]
    )
    assert replay_wls_interval_certificate(certificate, execution, digital)["status"] == "MATCH"


def test_dyadic_radius_is_exactly_the_published_float_and_is_positive() -> None:
    execution, digital = _execution()
    certificate = build_wls_interval_certificate(execution, digital, precision_bits=192)
    floors = certified_radius_floors(certificate, execution, digital)

    for row in certificate["response_disks"]:
        radius = floors[row["harmonic"]]
        assert radius == Fraction.from_float(row["absolute_error_upper"]["float_upper"])
        assert radius > 0
        assert radius < Fraction(1, 1_000_000)


def test_nonidentity_complex_front_end_is_included_as_exact_binary64() -> None:
    front_end = {1: 0.8 + 0.1j, 2: 0.7 - 0.2j, 3: 0.5 + 0.25j}
    execution, digital = _execution(front_end=front_end)

    certificate = build_wls_interval_certificate(execution, digital, precision_bits=192)

    assert certificate["relation"] == "CERTIFIED_INTERVAL"
    rows = certificate["mathematical_contract"]["front_end_exact_dyadics"]
    assert [row["harmonic"] for row in rows] == [1, 2, 3]
    assert Fraction(
        int(rows[0]["real"]["numerator"]), int(rows[0]["real"]["denominator"])
    ) == Fraction.from_float(front_end[1].real)
    assert replay_wls_interval_certificate(certificate, execution, digital)["status"] == "MATCH"


def test_issuance_requires_execution_match_and_exact_payload() -> None:
    execution, digital = _execution()
    changed = digital.copy()
    changed[5] += 1

    with pytest.raises(ValueError, match="execution replay must be MATCH"):
        build_wls_interval_certificate(execution, changed)
    with pytest.raises(ValueError, match="exact integer dtype"):
        build_wls_interval_certificate(execution, digital.astype(np.float64))
    with pytest.raises(ValueError, match="precision"):
        build_wls_interval_certificate(execution, digital, precision_bits=True)


def test_execution_numeric_subclass_is_rejected_before_dyadic_conversion() -> None:
    execution, digital = _execution()
    stored = execution["result"]["z_tilde"][0]["real"]
    execution["result"]["z_tilde"][0]["real"] = _RatioZero(stored)

    with pytest.raises(ValueError, match="exact plain tree"):
        build_wls_interval_certificate(execution, digital)


def test_payload_validation_takes_a_private_int64_snapshot() -> None:
    original = np.arange(32, dtype=np.int64)
    snapshot = _validated_int64_payload(original, expected_count=original.size)

    assert not np.shares_memory(snapshot, original)
    original[0] = 999
    assert snapshot[0] == 0


def test_replay_detects_certificate_and_execution_tampering() -> None:
    execution, digital = _execution()
    certificate = build_wls_interval_certificate(execution, digital, precision_bits=192)

    tampered = copy.deepcopy(certificate)
    tampered["response_disks"][0]["absolute_error_upper"]["float_upper"] *= 2.0
    assert replay_wls_interval_certificate(tampered, execution, digital) == {
        "status": "MISMATCH",
        "reason": "stored_certificate_hash_mismatch",
    }

    resigned_execution = copy.deepcopy(execution)
    resigned_execution["sufficient_statistics"]["b"][0]["real"] += 0.01
    unsigned_execution = dict(resigned_execution)
    unsigned_execution.pop("execution_sha256")
    resigned_execution["execution_sha256"] = canonical_hash(unsigned_execution)
    replay = replay_wls_interval_certificate(certificate, resigned_execution, digital)
    assert replay["status"] == "UNKNOWN"


def test_replay_is_resource_bounded_and_global_precision_is_restored() -> None:
    execution, digital = _execution()
    previous = ctx.prec
    certificate = build_wls_interval_certificate(execution, digital, precision_bits=192)
    assert ctx.prec == previous
    replay_wls_interval_certificate(certificate, execution, digital)
    assert ctx.prec == previous

    deep = copy.deepcopy(certificate)
    cursor = deep
    for _ in range(40):
        cursor["nested"] = {}
        cursor = cursor["nested"]
    assert replay_wls_interval_certificate(deep, execution, digital) == {
        "status": "UNKNOWN",
        "reason": "document_depth_budget_exceeded",
    }

    numeric_subclass = copy.deepcopy(certificate)
    numeric_subclass["numeric_proof"]["precision_bits"] = np.int64(192)
    assert replay_wls_interval_certificate(numeric_subclass, execution, digital) == {
        "status": "UNKNOWN",
        "reason": "document_non_plain_json_value",
    }


def test_claim_scope_excludes_non_numeric_ecg_errors() -> None:
    execution, digital = _execution()
    certificate = build_wls_interval_certificate(execution, digital, precision_bits=192)

    assert certificate["claim_scope"] == "frozen_mathematical_sample_functional_only"
    assert set(certificate["exclusions"]) >= {
        "acquisition_and_sensor_error",
        "beat_annotation_and_window_timing_error",
        "anti_alias_error",
        "front_end_calibration_uncertainty",
        "model_discrepancy_and_physiological_interpretation",
    }


def test_zero_signal_radius_does_not_underflow_to_an_unproved_zero() -> None:
    execution, digital = _execution()
    zeros = np.zeros_like(digital)
    zero_execution = build_wls_execution_record(
        zeros,
        tuple(execution["window"]["beat_samples"]),
        sampling_frequency=250.0,
        adc_gain=1000.0,
        baseline=0,
        harmonics=(1, 2, 3),
        protocol_sha256="1" * 64,
        record_sha256="2" * 64,
        window_sha256="3" * 64,
    )

    certificate = build_wls_interval_certificate(zero_execution, zeros, precision_bits=192)
    floors = certified_radius_floors(certificate, zero_execution, zeros)

    assert certificate["relation"] == "CERTIFIED_INTERVAL"
    assert all(radius == Fraction.from_float(5e-324) for radius in floors.values())


def test_exactly_aliased_gram_fails_closed_even_if_binary64_execution_solved() -> None:
    beats = tuple(map(int, np.linspace(0, 12, 9, dtype=int)))
    digital = np.arange(12, dtype=np.int64)
    execution = build_wls_execution_record(
        digital,
        beats,
        sampling_frequency=1.0,
        adc_gain=1.0,
        baseline=0,
        harmonics=(1, 2),
        protocol_sha256="1" * 64,
        record_sha256="2" * 64,
        window_sha256="3" * 64,
    )

    certificate = build_wls_interval_certificate(execution, digital, precision_bits=192)

    assert certificate["relation"] == "UNKNOWN"
    assert certificate["reason"] == "interval_computation_undecided:ArithmeticError"
    assert certificate["response_disks"] == []
    assert replay_wls_interval_certificate(certificate, execution, digital) == {
        "status": "MATCH",
        "relation": "UNKNOWN",
        "reason": "interval_computation_undecided:ArithmeticError",
        "certificate_sha256": certificate["certificate_sha256"],
    }
    with pytest.raises(ValueError, match="replayed CERTIFIED_INTERVAL"):
        certified_radius_floors(certificate, execution, digital)
