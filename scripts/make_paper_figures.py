"""Generate deterministic figures from the locked Fantasia pilot bundle."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def make_fantasia_summary(bundle_dir: Path, output_dir: Path) -> Path:
    """Create the two-panel pilot figure used by the submission draft."""

    summary = json.loads((bundle_dir / "summary.json").read_text(encoding="utf-8"))
    model_fits = json.loads((bundle_dir / "model_fits.json").read_text(encoding="utf-8"))
    comparisons = summary["validation_comparisons"]
    record_order = ("f1y01", "f1o01", "f2y01", "f2o01")
    alpha_by_record = {
        row["record"]: row["models"]["fractional_common"]["alpha"] for row in model_fits["records"]
    }

    constant_delta = np.asarray(
        [
            100.0 * (comparisons[record]["fractional_vs_constant_ratio"] - 1.0)
            for record in record_order
        ]
    )
    shuffle_delta = np.asarray(
        [
            100.0 * (comparisons[record]["fractional_vs_rate_shuffle_ratio"] - 1.0)
            for record in record_order
        ]
    )
    alpha = np.asarray([alpha_by_record[record] for record in record_order])

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    figure, axes = plt.subplots(1, 2, figsize=(7.2, 2.8), constrained_layout=True)
    x = np.arange(len(record_order))
    width = 0.34
    axes[0].axhline(0.0, color="0.25", linewidth=0.8)
    axes[0].bar(
        x - width / 2,
        constant_delta,
        width,
        label="vs constant Fourier",
        color="#4472C4",
    )
    axes[0].bar(
        x + width / 2,
        shuffle_delta,
        width,
        label="vs rate shuffle",
        color="#ED7D31",
    )
    axes[0].set_xticks(x, record_order)
    axes[0].set_ylabel("Fractional RMSE difference (%)")
    axes[0].set_title("(a) Held-out comparison")
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    marker_colors = ["#A5A5A5", "#A5A5A5", "#70AD47", "#70AD47"]
    axes[1].axhline(1.0, color="#C00000", linewidth=1.0)
    axes[1].scatter(x, alpha, s=48, c=marker_colors, edgecolor="black", linewidth=0.5, zorder=3)
    axes[1].set_xticks(x, record_order)
    axes[1].set_ylim(0.94, 1.01)
    axes[1].set_ylabel(r"Fitted fractional order $\alpha$")
    axes[1].set_title("(b) Boundary diagnostic")
    axes[1].text(
        0.04,
        0.14,
        "red line: $\\alpha=1$ boundary\n"
        "4/4 returned solutions at the nested boundary\n"
        "(not a global-optimum certificate)",
        transform=axes[1].transAxes,
        fontsize=7.5,
        va="bottom",
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / "fantasia_pilot_summary.pdf"
    figure.savefig(
        destination,
        bbox_inches="tight",
        metadata={
            "Creator": "FO-EKF locked figure generator",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    plt.close(figure)
    return destination


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bundle-dir",
        type=Path,
        default=root / "research" / "results" / "fantasia_pilot_2026-07-12",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "paper" / "figures",
    )
    arguments = parser.parse_args()
    destination = make_fantasia_summary(arguments.bundle_dir, arguments.output_dir)
    print(destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
