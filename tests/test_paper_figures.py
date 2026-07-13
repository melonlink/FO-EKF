import shutil
from pathlib import Path

import pytest

from scripts.make_paper_figures import (
    make_certification_workflow_comparison,
    make_core_method_evidence,
    make_fantasia_summary,
    make_identifiability_geometry,
)


def test_locked_fantasia_figure_is_generated(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    bundle = root / "research" / "results" / "fantasia_pilot_2026-07-12"
    destination = make_fantasia_summary(
        bundle,
        tmp_path,
    )

    assert destination == tmp_path / "fantasia_pilot_summary.pdf"
    assert destination.read_bytes().startswith(b"%PDF-")
    assert destination.stat().st_size > 5_000

    tampered_bundle = tmp_path / "tampered_fantasia"
    shutil.copytree(bundle, tampered_bundle)
    with (tampered_bundle / "model_fits.json").open("ab") as stream:
        stream.write(b"\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        make_fantasia_summary(tampered_bundle, tmp_path)


def test_method_figures_are_generated_from_the_locked_evidence(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    destinations = (
        make_identifiability_geometry(tmp_path),
        make_certification_workflow_comparison(tmp_path),
        make_core_method_evidence(
            root / "research" / "results" / "two_rate_recovery_2026-07-12",
            root / "research" / "results" / "synthetic_certificate_chain_2026-07-12",
            root / "research" / "results" / "s0_s3_2026-07-11",
            tmp_path,
        ),
    )

    assert {path.name for path in destinations} == {
        "method_identifiability_geometry.pdf",
        "certification_workflow_comparison.pdf",
        "core_method_evidence.pdf",
    }
    for destination in destinations:
        assert destination.read_bytes().startswith(b"%PDF-")
        assert destination.stat().st_size > 5_000

    tampered_two_rate = tmp_path / "tampered_two_rate"
    shutil.copytree(
        root / "research" / "results" / "two_rate_recovery_2026-07-12",
        tampered_two_rate,
    )
    with (tampered_two_rate / "cases.csv").open("ab") as stream:
        stream.write(b"\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        make_core_method_evidence(
            tampered_two_rate,
            root / "research" / "results" / "synthetic_certificate_chain_2026-07-12",
            root / "research" / "results" / "s0_s3_2026-07-11",
            tmp_path,
        )

    tampered_synthetic = tmp_path / "tampered_synthetic"
    shutil.copytree(
        root / "research" / "results" / "synthetic_certificate_chain_2026-07-12",
        tampered_synthetic,
    )
    with (tampered_synthetic / "summary.json").open("ab") as stream:
        stream.write(b"\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        make_core_method_evidence(
            root / "research" / "results" / "two_rate_recovery_2026-07-12",
            tampered_synthetic,
            root / "research" / "results" / "s0_s3_2026-07-11",
            tmp_path,
        )
