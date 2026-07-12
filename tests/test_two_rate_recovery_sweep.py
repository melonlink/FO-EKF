import hashlib
import json
from pathlib import Path

import pytest

import scripts.run_two_rate_recovery_sweep as sweep

FROZEN_COMMIT = "a" * 40
RESULT_COMMIT = "b" * 40


def _generate_test_evidence(output_dir: Path) -> dict[str, object]:
    return sweep.run_two_rate_recovery_sweep(output_dir, git_commit=FROZEN_COMMIT)


def test_deterministic_two_rate_and_three_rate_sweep_reuses_frozen_commit(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    summary = _generate_test_evidence(tmp_path)
    first_bytes = {
        name: (tmp_path / name).read_bytes() for name in ("cases.csv", "summary.json", "SHA256SUMS")
    }
    monkeypatch.setattr(sweep, "_production_sources_dirty", lambda: False)
    monkeypatch.setattr(sweep, "_current_git_commit", lambda: RESULT_COMMIT)
    repeated = sweep.run_two_rate_recovery_sweep(tmp_path)

    assert summary["design"]["case_count"] == 1260
    assert summary["design"]["seed"] is None
    assert summary["design"]["order_values"][0] == 0.05
    assert summary["design"]["order_values"][-1] == 1.0
    assert summary["design"]["damping_values"] == [0.01, 0.1, 1.0, 10.0, 100.0]
    assert summary["claim_boundary"]["analytical_proof_replaced"] is False
    assert summary["determinism"] == {
        "case_serialization_repeats": 2,
        "case_serialization_byte_identical": True,
        "fixed_grid_no_prng": True,
    }
    assert len(summary["provenance"]["production_source_manifest_sha256"]) == 64
    assert summary["provenance"]["environment"]["float_mantissa_bits"] == 53
    assert summary["provenance"]["code_commit"] == FROZEN_COMMIT
    assert repeated["provenance"]["code_commit"] == FROZEN_COMMIT
    assert any(
        entry["path"] == "requirements-repro-lock.txt"
        for entry in summary["provenance"]["production_source_manifest"]
    )

    two_rate = summary["outcomes"]["two_rate_global_inverse"]
    assert two_rate["successful_cases"] == 1260
    assert two_rate["failed_cases"] == 0
    assert two_rate["maximum_order_absolute_error"] < 1.0e-8
    assert two_rate["maximum_damping_absolute_error"] < 1.0e-7

    invariance = summary["outcomes"]["common_complex_gain_invariance"]
    assert invariance["successful_cases"] == 1260
    assert invariance["failed_cases"] == 0
    assert invariance["maximum_order_difference"] < 1.0e-8
    assert invariance["maximum_damping_difference"] < 1.0e-7

    positive = summary["outcomes"]["three_rate_positive_consistency"]
    assert positive["consistent_cases"] == 1260
    assert positive["inconsistent_cases"] == 0
    assert positive["maximum_geometric_cross_ratio_absolute_error"] < 1.0e-10

    tamper = summary["outcomes"]["deterministic_third_rate_tamper"]
    assert tamper["detected_cases"] == 1260
    assert tamper["missed_cases"] == 0
    assert tamper["consistency_outcome_counts"] == {"INCONSISTENT": 1260}
    assert tamper["minimum_tampered_cross_ratio_deviation"] > 1.0e-8
    assert len((tmp_path / "cases.csv").read_text(encoding="utf-8").splitlines()) == 1261

    assert repeated["summary_sha256"] == summary["summary_sha256"]
    for name, contents in first_bytes.items():
        assert (tmp_path / name).read_bytes() == contents

    cases_sha256 = hashlib.sha256((tmp_path / "cases.csv").read_bytes()).hexdigest()
    summary_sha256 = hashlib.sha256((tmp_path / "summary.json").read_bytes()).hexdigest()
    assert summary["artifact_files"]["cases.csv"] == cases_sha256
    assert (tmp_path / "SHA256SUMS").read_text(encoding="ascii").splitlines() == [
        f"{cases_sha256}  cases.csv",
        f"{summary_sha256}  summary.json",
    ]

    unsealed = json.loads((tmp_path / "summary.json").read_text(encoding="utf-8"))
    stored_seal = unsealed.pop("summary_sha256")
    canonical = json.dumps(
        unsealed,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    assert stored_seal == hashlib.sha256(canonical).hexdigest()


def test_formal_generation_refuses_dirty_production_sources(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(sweep, "_production_sources_dirty", lambda: True)

    with pytest.raises(RuntimeError, match="dirty production sources"):
        sweep.run_two_rate_recovery_sweep(tmp_path)

    assert not tmp_path.exists() or not any(tmp_path.iterdir())


def test_existing_evidence_rejects_changed_production_manifest(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _generate_test_evidence(tmp_path)
    monkeypatch.setattr(sweep, "_production_sources_dirty", lambda: False)
    monkeypatch.setattr(sweep, "_current_git_commit", lambda: RESULT_COMMIT)
    monkeypatch.setattr(
        sweep,
        "_production_source_manifest",
        lambda: (
            "0" * 64,
            [{"path": "changed-production-source", "sha256": "0" * 64}],
        ),
    )

    with pytest.raises(RuntimeError, match="different production sources"):
        sweep.run_two_rate_recovery_sweep(tmp_path)


@pytest.mark.parametrize(
    ("target", "expected_error"),
    (
        ("summary.json", "invalid internal seal"),
        ("cases.csv", "case artifact hash mismatch"),
        ("SHA256SUMS", "invalid or duplicate entry"),
    ),
)
def test_existing_evidence_tamper_is_rejected(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    target: str,
    expected_error: str,
) -> None:
    _generate_test_evidence(tmp_path)
    monkeypatch.setattr(sweep, "_production_sources_dirty", lambda: False)
    monkeypatch.setattr(sweep, "_current_git_commit", lambda: RESULT_COMMIT)
    path = tmp_path / target
    if target == "summary.json":
        document = json.loads(path.read_text(encoding="utf-8"))
        document["design"]["case_count"] = 1259
        path.write_text(json.dumps(document), encoding="utf-8")
    elif target == "cases.csv":
        path.write_bytes(path.read_bytes() + b"tamper\n")
    else:
        first_line = path.read_text(encoding="ascii").splitlines()[0]
        path.write_text(f"{first_line}\n{first_line}\n", encoding="ascii")

    with pytest.raises(RuntimeError, match=expected_error):
        sweep.run_two_rate_recovery_sweep(tmp_path)
