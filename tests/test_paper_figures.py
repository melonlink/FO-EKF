from pathlib import Path

from scripts.make_paper_figures import make_fantasia_summary


def test_locked_fantasia_figure_is_generated(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    destination = make_fantasia_summary(
        root / "research" / "results" / "fantasia_pilot_2026-07-12",
        tmp_path,
    )

    assert destination == tmp_path / "fantasia_pilot_summary.pdf"
    assert destination.read_bytes().startswith(b"%PDF-")
    assert destination.stat().st_size > 5_000
