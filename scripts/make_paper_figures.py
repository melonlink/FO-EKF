"""Generate deterministic vector figures for the FO-ECG submission."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

BLUE = "#0072B2"
SKY = "#56B4E9"
GREEN = "#009E73"
ORANGE = "#E69F00"
VERMILLION = "#D55E00"
PURPLE = "#CC79A7"
DARK = "#30343B"
LIGHT_GRAY = "#F2F3F5"
MID_GRAY = "#7A7F87"


def _paper_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def _save_vector(figure: plt.Figure, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
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


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _verify_sealed_bundle(directory: Path, expected_files: set[str]) -> dict[str, str]:
    """Fail closed unless the tracked SHA256SUMS exactly seals the expected files."""

    manifest_path = directory / "SHA256SUMS"
    entries: dict[str, str] = {}
    for line in manifest_path.read_text(encoding="ascii").splitlines():
        parts = line.split("  ", maxsplit=1)
        if len(parts) != 2 or len(parts[0]) != 64 or parts[1] in entries:
            raise ValueError(f"malformed sealed bundle manifest: {manifest_path}")
        entries[parts[1]] = parts[0]
    if set(entries) != expected_files:
        raise ValueError(f"sealed bundle allowlist mismatch: {directory}")
    for name, expected in entries.items():
        if _sha256_file(directory / name) != expected:
            raise ValueError(f"sealed bundle checksum mismatch: {directory / name}")
    return entries


def make_fantasia_summary(bundle_dir: Path, output_dir: Path) -> Path:
    """Create the two-panel pilot figure used by the submission draft."""

    _verify_sealed_bundle(
        bundle_dir,
        {
            "design_lock.json",
            "heldout_metrics.csv",
            "model_fits.json",
            "r2_status.json",
            "source_inventory.json",
            "summary.json",
            "window_manifest.csv.gz",
            "wls_coefficients.csv.gz",
            "wls_execution_records.jsonl.gz",
        },
    )
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

    _paper_style()
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
    figure.text(
        0.5,
        -0.035,
        "Claim gate: fractional specificity NOT SUPPORTED   |   "
        "Coverage gate: real-data R2 NOT CERTIFIABLE",
        ha="center",
        va="top",
        fontsize=8.0,
        color=VERMILLION,
        weight="bold",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "#FBE7DE", "edgecolor": VERMILLION},
    )

    destination = output_dir / "fantasia_pilot_summary.pdf"
    return _save_vector(figure, destination)


def make_identifiability_geometry(output_dir: Path) -> Path:
    """Draw the one-/two-/three-rate identifiability geometry."""

    _paper_style()
    figure, axes = plt.subplots(1, 3, figsize=(7.25, 2.55), constrained_layout=True)

    # One rate: q can absorb every admissible (alpha, lambda) pair.
    axis = axes[0]
    axis.set_facecolor("#F7F7F7")
    alpha_grid = np.linspace(0.08, 1.0, 18)
    damping_grid = np.linspace(-2.0, 2.0, 15)
    aa, ll = np.meshgrid(alpha_grid, damping_grid)
    axis.scatter(aa, ll, s=5, color="#C7CBD1", alpha=0.65, linewidths=0)
    feasible = np.asarray([(0.16, -1.4), (0.32, 1.1), (0.55, -0.2), (0.76, 1.65), (0.94, 0.45)])
    axis.scatter(
        feasible[:, 0],
        feasible[:, 1],
        s=31,
        c=[BLUE, ORANGE, GREEN, PURPLE, VERMILLION],
        edgecolor="white",
        linewidth=0.6,
        zorder=3,
    )
    axis.set_xlim(0.04, 1.03)
    axis.set_ylim(-2.15, 2.15)
    axis.set_xlabel(r"order $\alpha$")
    axis.set_ylabel(r"$\log_{10}\lambda$")
    axis.set_title("(a) One rate\nnon-identifiable", fontsize=8.8, pad=6)
    axis.text(
        0.04,
        0.97,
        r"$q(\alpha,\lambda)=Z_\star[\lambda+(i\nu)^\alpha]$" + "\n"
        r"every dot reproduces the same $Z_\star$",
        transform=axis.transAxes,
        va="top",
        fontsize=7.1,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#D0D3D8"},
    )

    # Two rates: quotient rays are independent of the common complex gauge g*q.
    axis = axes[1]
    rate_ratio = 2.0
    orders = np.asarray([0.30, 0.45, 0.60, 0.75, 0.90])
    colors = plt.colormaps["viridis"](np.linspace(0.12, 0.88, len(orders)))
    damping_ratios = np.linspace(0.02, 1.35, 100)
    for order, color in zip(orders, colors, strict=True):
        denominator = rate_ratio**order - 1.0
        direction = np.exp(-0.5j * np.pi * order)
        quotient = (1.0 + damping_ratios * direction) / denominator
        axis.plot(quotient.real, quotient.imag, color=color, linewidth=1.25)
        axis.text(
            quotient.real[-1],
            quotient.imag[-1],
            rf" ${order:.2f}$",
            color=color,
            fontsize=6.1,
            va="center",
        )
    observed_order = 0.60
    observed_damping_ratio = 0.85
    observed_denominator = rate_ratio**observed_order - 1.0
    observed = (
        1.0 + observed_damping_ratio * np.exp(-0.5j * np.pi * observed_order)
    ) / observed_denominator
    axis.plot(
        observed.real,
        observed.imag,
        marker="*",
        markersize=9,
        markerfacecolor=VERMILLION,
        markeredgecolor="black",
        markeredgewidth=0.55,
        zorder=5,
    )
    axis.annotate(
        "unique ray",
        xy=(observed.real, observed.imag),
        xytext=(4.4, -0.75),
        fontsize=6.8,
        arrowprops={"arrowstyle": "->", "color": DARK, "linewidth": 0.8},
    )
    axis.axhline(0.0, color="#B8BCC2", linewidth=0.7)
    axis.set_xlim(0.75, 10.2)
    axis.set_ylim(-4.2, 0.25)
    axis.set_xlabel(r"$\operatorname{Re} C_{12}$")
    axis.set_ylabel(r"$\operatorname{Im} C_{12}$")
    axis.set_title("(b) Two rates\nexact global inverse", fontsize=8.8, pad=6)
    axis.text(
        0.04,
        0.97,
        r"$W_r=1/Z_r,\ C_{12}=W_1/(W_2-W_1)$" + "\n"
        r"exact/noiseless; common complex $gq$ cancels",
        transform=axis.transAxes,
        va="top",
        fontsize=7.1,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#D0D3D8"},
    )
    axis.text(
        0.91,
        0.08,
        r"ray label: $\alpha$",
        transform=axis.transAxes,
        ha="right",
        fontsize=6.2,
    )

    # Three geometric rates: the model predicts a real cross-ratio segment.
    axis = axes[2]
    axis.axhline(0.0, color="#B8BCC2", linewidth=0.7)
    axis.plot([2.0, 3.0], [0.0, 0.0], color=GREEN, linewidth=5.0, solid_capstyle="round")
    axis.plot(2.0, 0.0, marker="o", markersize=5.5, markerfacecolor="white", markeredgecolor=GREEN)
    axis.plot(3.0, 0.0, marker="o", markersize=5.5, markerfacecolor=GREEN, markeredgecolor=GREEN)
    valid_cross_ratio = 1.0 + rate_ratio**observed_order
    axis.plot(
        valid_cross_ratio,
        0.0,
        marker="o",
        markersize=6.5,
        markerfacecolor=BLUE,
        markeredgecolor="black",
        markeredgewidth=0.5,
        zorder=4,
    )
    tampered = valid_cross_ratio + 0.16 + 0.28j
    axis.plot(
        tampered.real,
        tampered.imag,
        marker="X",
        markersize=7,
        color=VERMILLION,
        markeredgecolor="black",
        markeredgewidth=0.4,
        zorder=4,
    )
    axis.annotate(
        "off-manifold third rate\nfalsifies exact model",
        xy=(tampered.real, tampered.imag),
        xytext=(2.74, 0.12),
        fontsize=6.6,
        ha="center",
        arrowprops={"arrowstyle": "->", "color": DARK, "linewidth": 0.8},
    )
    axis.text(valid_cross_ratio, -0.12, "consistent", color=BLUE, fontsize=6.6, ha="center")
    axis.set_xlim(1.82, 3.24)
    axis.set_ylim(-0.48, 0.58)
    axis.set_xlabel(r"$\operatorname{Re} Q_{123}$")
    axis.set_ylabel(r"$\operatorname{Im} Q_{123}$")
    axis.set_title("(c) Three rates\nexact consistency test", fontsize=8.8, pad=6)
    axis.text(
        0.04,
        0.97,
        r"$Q_{123}=\dfrac{W_3-W_1}{W_2-W_1}=1+r^\alpha\in(2,3]$" + "\n"
        r"geometric rates, $r=2$; disks need a joint-set test",
        transform=axis.transAxes,
        va="top",
        fontsize=6.9,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#D0D3D8"},
    )

    destination = output_dir / "method_identifiability_geometry.pdf"
    return _save_vector(figure, destination)


def _workflow_box(
    axis: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    text: str,
    *,
    facecolor: str,
    edgecolor: str,
    fontsize: float = 7.0,
) -> None:
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.008,rounding_size=0.012",
        linewidth=1.0,
        facecolor=facecolor,
        edgecolor=edgecolor,
    )
    axis.add_patch(patch)
    axis.text(x + width / 2.0, y + height / 2.0, text, ha="center", va="center", fontsize=fontsize)


def _workflow_arrow(axis: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    axis.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=9,
            linewidth=1.0,
            color=MID_GRAY,
            shrinkA=2,
            shrinkB=2,
        )
    )


def make_certification_workflow_comparison(output_dir: Path) -> Path:
    """Contrast a point-estimation path with the fail-closed certificate path."""

    _paper_style()
    figure, axis = plt.subplots(figsize=(7.25, 3.15), constrained_layout=True)
    axis.set_xlim(0.0, 1.0)
    axis.set_ylim(0.0, 1.0)
    axis.axis("off")

    axis.add_patch(Rectangle((0.0, 0.58), 1.0, 0.36, facecolor="#F7F7F7", edgecolor="none"))
    axis.add_patch(Rectangle((0.0, 0.20), 1.0, 0.34, facecolor="#F1F7FA", edgecolor="none"))
    axis.text(
        0.012,
        0.895,
        "Conventional point inversion (illustrative; not evaluated here)",
        color=DARK,
        weight="bold",
        fontsize=8.3,
    )
    axis.text(0.012, 0.505, "Certified set path", color=BLUE, weight="bold", fontsize=8.7)

    top_x = (0.13, 0.355, 0.58, 0.805)
    top_y = 0.675
    top_width = 0.165
    top_height = 0.145
    top_labels = (
        "Integer samples\n$y_n$",
        "Point WLS\n$\\hat Z_k$",
        "Point optimizer / filter\npoint state $\\hat\\theta$",
        "Point output\nresidual / status",
    )
    top_faces = ("white", "#E8F3F8", "#EFEAF6", "#EEEEEE")
    top_edges = (MID_GRAY, BLUE, PURPLE, MID_GRAY)
    for x, label, face, edge in zip(top_x, top_labels, top_faces, top_edges, strict=True):
        _workflow_box(
            axis,
            x,
            top_y,
            top_width,
            top_height,
            label,
            facecolor=face,
            edgecolor=edge,
        )
    for left, right in zip(top_x, top_x[1:], strict=False):
        _workflow_arrow(
            axis,
            (left + top_width, top_y + top_height / 2.0),
            (right, top_y + top_height / 2.0),
        )
    axis.text(
        0.5,
        0.615,
        "Useful point estimate; not an end-to-end enclosure of acquisition, WLS, and model error.",
        ha="center",
        fontsize=6.8,
        color="#555A61",
    )

    input_x, input_y, input_width, input_height = 0.015, 0.315, 0.12, 0.145
    wls_x, wls_y, wls_width, wls_height = 0.165, 0.397, 0.16, 0.105
    primitive_x, primitive_y = 0.165, 0.245
    gate_x, gate_y, gate_width, gate_height = 0.365, 0.315, 0.13, 0.145
    response_x, response_y, response_width, response_height = 0.535, 0.315, 0.125, 0.145
    joint_x, joint_y, joint_width, joint_height = 0.695, 0.315, 0.125, 0.145
    output_x, output_y, output_width, output_height = 0.855, 0.315, 0.13, 0.145
    _workflow_box(
        axis,
        input_x,
        input_y,
        input_width,
        input_height,
        "Signed int64\n+ exact indices\n+ frozen hashes",
        facecolor="white",
        edgecolor=DARK,
        fontsize=6.5,
    )
    _workflow_box(
        axis,
        wls_x,
        wls_y,
        wls_width,
        wls_height,
        "Replayed WLS centre\n+ Arb/ACB numerical disk",
        facecolor="#DDEFF7",
        edgecolor=BLUE,
        fontsize=6.5,
    )
    _workflow_box(
        axis,
        primitive_x,
        primitive_y,
        wls_width,
        wls_height,
        "Independent R2 bounds\n$\\rho_{\\rm ADC},\\rho_t,\\rho_{\\rm model},\\ldots$",
        facecolor="#FFF1CF",
        edgecolor=ORANGE,
        fontsize=6.2,
    )
    _workflow_box(
        axis,
        gate_x,
        gate_y,
        gate_width,
        gate_height,
        "R2 coverage gate\nhashes + bounds + Gram",
        facecolor="#FFF7E5",
        edgecolor=ORANGE,
        fontsize=6.2,
    )
    _workflow_box(
        axis,
        response_x,
        response_y,
        response_width,
        response_height,
        "R2 response disks\n$Z_k\\in\\mathbb{D}_k$",
        facecolor="#E8F3F8",
        edgecolor=SKY,
        fontsize=6.7,
    )
    _workflow_box(
        axis,
        joint_x,
        joint_y,
        joint_width,
        joint_height,
        "Joint inversion\nshared $\\alpha,\\lambda,q_m$",
        facecolor="#E1F2EC",
        edgecolor=GREEN,
        fontsize=6.7,
    )
    _workflow_box(
        axis,
        output_x,
        output_y,
        output_width,
        output_height,
        "Certified order set\n/ preserved witness",
        facecolor="#DDF2E9",
        edgecolor=GREEN,
        fontsize=6.1,
    )
    _workflow_arrow(
        axis,
        (input_x + input_width, input_y + input_height / 2.0),
        (wls_x, wls_y + wls_height / 2.0),
    )
    _workflow_arrow(
        axis,
        (wls_x + wls_width, wls_y + wls_height / 2.0),
        (gate_x, gate_y + 0.68 * gate_height),
    )
    _workflow_arrow(
        axis,
        (primitive_x + wls_width, primitive_y + wls_height / 2.0),
        (gate_x, gate_y + 0.32 * gate_height),
    )
    _workflow_arrow(
        axis,
        (gate_x + gate_width, gate_y + gate_height / 2.0),
        (response_x, response_y + response_height / 2.0),
    )
    axis.text(
        (gate_x + gate_width + response_x) / 2.0,
        0.468,
        r"$\mathrm{PASS}_{*}$",
        ha="center",
        fontsize=6.0,
        color=GREEN,
    )
    _workflow_arrow(
        axis,
        (response_x + response_width, response_y + response_height / 2.0),
        (joint_x, joint_y + joint_height / 2.0),
    )
    _workflow_arrow(
        axis,
        (joint_x + joint_width, joint_y + joint_height / 2.0),
        (output_x, output_y + output_height / 2.0),
    )

    nc_x, nc_y, nc_width, nc_height = 0.345, 0.14, 0.17, 0.085
    _workflow_box(
        axis,
        nc_x,
        nc_y,
        nc_width,
        nc_height,
        "NOT CERTIFIABLE / EXCLUDE\nmissing bound / invalid protocol",
        facecolor="#FBE7DE",
        edgecolor=VERMILLION,
        fontsize=5.8,
    )
    _workflow_arrow(
        axis,
        (gate_x + gate_width / 2.0, gate_y),
        (nc_x + nc_width / 2.0, nc_y + nc_height),
    )

    axis.text(
        0.735,
        0.275,
        "R2 PASS forms disks; joint inversion preserves shared parameters across rates.",
        ha="center",
        fontsize=6.2,
        color=BLUE,
        weight="bold",
    )
    axis.text(
        0.5,
        0.065,
        "R2 PASS means coverage, not model acceptance. REJECT MODEL requires exhaustive "
        "outer-set emptiness; no real-ECG PASS is asserted here.",
        ha="center",
        va="center",
        fontsize=6.5,
        color=DARK,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": LIGHT_GRAY, "edgecolor": "#CDD1D6"},
    )

    destination = output_dir / "certification_workflow_comparison.pdf"
    return _save_vector(figure, destination)


def make_core_method_evidence(
    two_rate_dir: Path,
    synthetic_dir: Path,
    legacy_dir: Path,
    output_dir: Path,
) -> Path:
    """Plot evidence aligned with the inverse, falsification, and disk claims."""

    two_rate_seal = _verify_sealed_bundle(two_rate_dir, {"cases.csv", "summary.json"})
    synthetic_seal = _verify_sealed_bundle(
        synthetic_dir,
        {"artifacts.json", "summary.json"},
    )
    cases = np.genfromtxt(
        two_rate_dir / "cases.csv",
        delimiter=",",
        names=True,
        dtype=None,
        encoding="utf-8",
    )
    two_rate_summary = json.loads((two_rate_dir / "summary.json").read_text(encoding="utf-8"))
    synthetic = json.loads((synthetic_dir / "summary.json").read_text(encoding="utf-8"))
    legacy = np.genfromtxt(
        legacy_dir / "s3_bound_summary.csv",
        delimiter=",",
        names=True,
        dtype=None,
        encoding="utf-8",
    )

    if two_rate_summary["artifact_files"]["cases.csv"] != two_rate_seal["cases.csv"]:
        raise ValueError("two-rate summary does not bind the sealed cases.csv")
    if synthetic["artifact_files"]["artifacts.json"] != synthetic_seal["artifacts.json"]:
        raise ValueError("synthetic summary does not bind the sealed artifacts.json")

    outcomes = two_rate_summary["outcomes"]
    recovery = outcomes["two_rate_global_inverse"]
    positive = outcomes["three_rate_positive_consistency"]
    tamper = outcomes["deterministic_third_rate_tamper"]
    case_count = len(cases)
    if recovery["successful_cases"] + recovery["failed_cases"] != case_count:
        raise ValueError("two-rate case count does not match the sealed summary")
    if positive["consistent_cases"] + positive["inconsistent_cases"] != case_count:
        raise ValueError("positive three-rate count does not match the sealed summary")
    if tamper["detected_cases"] + tamper["missed_cases"] != case_count:
        raise ValueError("tamper count does not match the sealed summary")
    attacks = synthetic["tamper_controls"]
    if len(attacks) != 2 or any(
        attack["relation"] != "NOT_CERTIFIABLE" for attack in attacks.values()
    ):
        raise ValueError("synthetic tamper controls are not the sealed fail-closed pair")

    _paper_style()
    figure, axes_grid = plt.subplots(2, 2, figsize=(7.25, 5.05), constrained_layout=True)
    axes = axes_grid.ravel()

    axis = axes[0]
    ratios = np.unique(cases["frequency_ratio"])
    ratio_colors = (BLUE, SKY, GREEN, PURPLE)
    for ratio, color in zip(ratios, ratio_colors, strict=True):
        selected = cases["frequency_ratio"] == ratio
        raw_errors = np.maximum(cases["raw_order_abs_error"][selected], 1.0e-16)
        gained_errors = np.maximum(cases["gained_order_abs_error"][selected], 1.0e-16)
        axis.scatter(
            cases["order_true"][selected],
            raw_errors,
            s=9,
            alpha=0.55,
            color=color,
            linewidths=0,
            label=rf"$r={ratio:g}$",
        )
        axis.scatter(
            cases["order_true"][selected],
            gained_errors,
            s=9,
            alpha=0.55,
            color=color,
            marker="x",
            linewidths=0.55,
        )
    maximum_order_error = float(
        max(np.max(cases["raw_order_abs_error"]), np.max(cases["gained_order_abs_error"]))
    )
    maximum_gain_difference = float(np.max(cases["common_gain_order_difference"]))
    axis.axhline(maximum_order_error, color=VERMILLION, linewidth=0.9, linestyle="--")
    axis.set_yscale("log")
    axis.set_ylim(8.0e-17, 1.2e-11)
    axis.set_xlim(0.03, 1.02)
    axis.set_xlabel(r"true order $\alpha$")
    axis.set_ylabel(r"absolute error $|\widehat\alpha-\alpha|$")
    axis.set_title("(a) Global inverse and common-gain audit", fontsize=9.2, pad=7)
    axis.legend(frameon=False, fontsize=6.2, ncol=2, loc="lower right")
    axis.text(
        0.04,
        0.96,
        f"{recovery['successful_cases']:,}/{case_count:,} original + gained\n"
        + rf"max. error $={maximum_order_error:.2e}$"
        + "\n"
        + rf"max. gain difference $={maximum_gain_difference:.2e}$",
        transform=axis.transAxes,
        va="top",
        fontsize=6.5,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#D0D3D8"},
    )
    axis.text(
        0.04,
        0.61,
        r"$\bullet$ original response; $\times$ after common complex gain"
        "\nplot floor $=10^{-16}$",
        transform=axis.transAxes,
        fontsize=6.1,
        color=DARK,
    )

    axis = axes[1]
    valid_residual = np.maximum(cases["geometric_cross_ratio_abs_error"], 1.0e-16)
    tamper_residual = np.maximum(cases["tampered_cross_ratio_deviation"], 1.0e-16)
    boxplot = axis.boxplot(
        [valid_residual, tamper_residual],
        tick_labels=("valid triple", "third-rate\ntamper"),
        patch_artist=True,
        showfliers=False,
        whis=(0, 100),
        widths=0.55,
        medianprops={"color": DARK, "linewidth": 1.0},
        whiskerprops={"color": MID_GRAY, "linewidth": 0.9},
        capprops={"color": MID_GRAY, "linewidth": 0.9},
    )
    for patch_box, color in zip(boxplot["boxes"], (GREEN, VERMILLION), strict=True):
        patch_box.set_facecolor(color)
        patch_box.set_alpha(0.45)
        patch_box.set_edgecolor(color)
    axis.set_yscale("log")
    axis.set_ylim(8.0e-17, max(1.0, 1.2 * float(np.max(tamper_residual))))
    axis.set_ylabel(r"$|Q_{\rm obs}-(1+r^{\alpha_\star})|$")
    axis.set_title("(b) Three-rate falsification", fontsize=9.2, pad=7)
    axis.text(
        0.5,
        0.96,
        f"{positive['consistent_cases']:,}/{case_count:,} consistent\n"
        f"{tamper['detected_cases']:,}/{case_count:,} detected",
        transform=axis.transAxes,
        ha="center",
        va="top",
        fontsize=6.8,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#D0D3D8"},
    )
    separation = float(np.min(tamper_residual) / np.max(valid_residual))
    axis.annotate(
        rf"${separation:.2e}\times$ separation",
        xy=(2.0, float(np.min(tamper_residual))),
        xytext=(1.55, 2.0e-7),
        ha="center",
        fontsize=6.4,
        color=DARK,
        arrowprops={"arrowstyle": "->", "color": DARK, "linewidth": 0.8},
    )
    axis.text(
        0.04,
        0.04,
        r"truth-referenced exact audit; plot floor $=10^{-16}$",
        transform=axis.transAxes,
        fontsize=6.0,
        color=DARK,
    )

    axis = axes[2]
    response = synthetic["response_disk"]
    radius = float(response["radius"])
    offset = (
        complex(
            float(response["known_model_response_real"]) - float(response["center_real"]),
            float(response["known_model_response_imag"]) - float(response["center_imag"]),
        )
        / radius
    )
    angle = np.linspace(0.0, 2.0 * np.pi, 400)
    axis.fill(np.cos(angle), np.sin(angle), color=SKY, alpha=0.22, zorder=1)
    axis.plot(np.cos(angle), np.sin(angle), color=BLUE, linewidth=1.3, zorder=2)
    axis.scatter(
        0.0,
        0.0,
        s=40,
        color=BLUE,
        edgecolor="black",
        linewidth=0.5,
        label="replayed centre",
        zorder=4,
    )
    axis.plot([0.0, offset.real], [0.0, offset.imag], color=ORANGE, linewidth=1.4, zorder=3)
    axis.scatter(
        offset.real,
        offset.imag,
        s=58,
        marker="*",
        color=GREEN,
        edgecolor="black",
        linewidth=0.5,
        label="known response",
        zorder=5,
    )
    normalized_error = abs(offset)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlim(-1.18, 1.18)
    axis.set_ylim(-1.18, 1.18)
    axis.set_xlabel(r"$(\operatorname{Re}Z-\operatorname{Re}\widetilde Z)/\epsilon$")
    axis.set_ylabel(r"$(\operatorname{Im}Z-\operatorname{Im}\widetilde Z)/\epsilon$")
    axis.set_title("(c) Exact-sample response disk", fontsize=9.2, pad=7)
    axis.legend(frameon=False, fontsize=6.1, loc="upper right")
    axis.text(
        0.03,
        0.04,
        rf"truth error $={normalized_error:.3f}\,\epsilon$" + "\n"
        rf"{len(attacks)} seal attacks $\rightarrow$ NOT CERTIFIABLE",
        transform=axis.transAxes,
        fontsize=6.7,
        va="bottom",
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": ORANGE},
    )

    axis = axes[3]
    relative_errors = np.unique(legacy["relative_error"])
    median_widths: list[float] = []
    target_pass_fractions: list[float] = []
    for relative_error in relative_errors:
        selected = legacy[legacy["relative_error"] == relative_error]
        median_widths.append(float(np.median(selected["median_chain_width"])))
        target_pass_fractions.append(
            float(np.sum(selected["primary_pass_0p02"]) / np.sum(selected["draws"]))
        )
    radius_percent = 100.0 * relative_errors
    axis.semilogx(
        radius_percent,
        median_widths,
        color=BLUE,
        marker="o",
        linewidth=1.6,
        label="median outer width",
    )
    axis.axhline(0.02, color=MID_GRAY, linestyle="--", linewidth=0.9, label="width target 0.02")
    axis.set_xlabel("relative response radius (%)")
    axis.set_ylabel(r"conservative order width $\Delta\alpha$")
    axis.set_ylim(0.0, 0.50)
    axis.set_title("(d) Finite-error resolution cost", fontsize=9.2, pad=7)
    second_axis = axis.twinx()
    second_axis.semilogx(
        radius_percent,
        target_pass_fractions,
        color=GREEN,
        marker="s",
        linewidth=1.4,
        label="fraction below target",
    )
    second_axis.set_ylabel(r"fraction with width $\leq0.02$", color=GREEN)
    second_axis.set_ylim(-0.02, 1.05)
    second_axis.tick_params(axis="y", colors=GREEN)
    handles, labels = axis.get_legend_handles_labels()
    handles_2, labels_2 = second_axis.get_legend_handles_labels()
    axis.legend(
        handles + handles_2,
        labels + labels_2,
        frameon=False,
        fontsize=6.1,
        loc="upper left",
    )
    axis.text(
        0.55,
        0.68,
        "legacy outer diagnostic\ncoverage retained",
        transform=axis.transAxes,
        ha="center",
        va="center",
        fontsize=6.1,
        color=DARK,
        bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": "#D0D3D8"},
    )

    destination = output_dir / "core_method_evidence.pdf"
    return _save_vector(figure, destination)


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
    parser.add_argument(
        "--two-rate-dir",
        type=Path,
        default=root / "research" / "results" / "two_rate_recovery_2026-07-12",
    )
    parser.add_argument(
        "--synthetic-dir",
        type=Path,
        default=root / "research" / "results" / "synthetic_certificate_chain_2026-07-12",
    )
    parser.add_argument(
        "--legacy-dir",
        type=Path,
        default=root / "research" / "results" / "s0_s3_2026-07-11",
    )
    arguments = parser.parse_args()
    destinations = (
        make_fantasia_summary(arguments.bundle_dir, arguments.output_dir),
        make_identifiability_geometry(arguments.output_dir),
        make_certification_workflow_comparison(arguments.output_dir),
        make_core_method_evidence(
            arguments.two_rate_dir,
            arguments.synthetic_dir,
            arguments.legacy_dir,
            arguments.output_dir,
        ),
    )
    for destination in destinations:
        print(destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
