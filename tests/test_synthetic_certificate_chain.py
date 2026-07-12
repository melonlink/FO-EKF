import hashlib
from pathlib import Path

from scripts.run_synthetic_certificate_chain import run_synthetic_chain


def test_synthetic_exact_sample_chain_closes_and_tampers_fail_closed(tmp_path: Path) -> None:
    summary = run_synthetic_chain(tmp_path, git_commit="a" * 40)
    first_summary_bytes = (tmp_path / "summary.json").read_bytes()
    repeated = run_synthetic_chain(tmp_path, git_commit="a" * 40)

    assert summary["stages"]["wls_execution"]["status"] == "MATCH"
    assert summary["stages"]["exact_interval"]["relation"] == "CERTIFIED_INTERVAL"
    assert summary["stages"]["r2_component"]["relation"] == "CERTIFIED_BOUND"
    assert summary["stages"]["protocol"]["status"] == "PASS_DETERMINISTIC"
    assert summary["stages"]["r2_t1_bundle"]["relation"] == "PRESERVED_ROBUST"
    assert summary["stages"]["r2_t1_bundle"]["replay_valid"] is True
    assert summary["numeric_radius_contract"]["solve_only_floor_added_again"] is False
    assert summary["coverage"]["mathematical_exact_sample_wls_coverage"] is True
    assert summary["coverage"]["uncovered_numeric_term"] is None
    assert summary["response_disk"]["disk_contains_known_response"] is True
    assert summary["response_disk"]["radius_to_response_magnitude_ratio"] < 0.01
    assert all(
        value["relation"] == "NOT_CERTIFIABLE" for value in summary["tamper_controls"].values()
    )
    assert (tmp_path / "summary.json").is_file()
    assert (tmp_path / "artifacts.json").is_file()
    assert (tmp_path / "SHA256SUMS").is_file()
    artifact_sha256 = hashlib.sha256((tmp_path / "artifacts.json").read_bytes()).hexdigest()
    summary_file_sha256 = hashlib.sha256((tmp_path / "summary.json").read_bytes()).hexdigest()
    assert summary["artifact_files"]["artifacts.json"] == artifact_sha256
    assert (tmp_path / "SHA256SUMS").read_text(encoding="ascii").splitlines() == [
        f"{artifact_sha256}  artifacts.json",
        f"{summary_file_sha256}  summary.json",
    ]
    assert repeated["summary_sha256"] == summary["summary_sha256"]
    assert (tmp_path / "summary.json").read_bytes() == first_summary_bytes
