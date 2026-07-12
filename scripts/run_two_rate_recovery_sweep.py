"""Run a deterministic numerical audit of the exact multi-rate inverse formulas.

This runner exercises production formulas from :mod:`fo_ekf.multirate` on a
fixed Cartesian grid.  It is a floating-point implementation audit of the
proved exact-model identities, not a replacement for their analytical proof
and not a noisy-ECG validation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from fo_ekf.multirate import (  # noqa: E402
    evaluate_reference_pair_consistency,
    harmonic_response,
    recover_reference_pair,
    response_invariant,
)

DEFAULT_OUTPUT = ROOT / "research" / "results" / "two_rate_recovery_2026-07-12"
SCHEMA = "fo-ekf.two-rate-recovery-sweep.v1"
PRODUCTION_SOURCE_PATHS = (
    ROOT / "pyproject.toml",
    ROOT / "requirements-repro-lock.txt",
    ROOT / "scripts" / "run_two_rate_recovery_sweep.py",
    ROOT / "src" / "fo_ekf" / "multirate.py",
)
CANONICAL_ARTIFACT_FILES = ("cases.csv", "summary.json", "SHA256SUMS")
ORDER_VALUES = tuple(0.05 + 0.95 * index / 20.0 for index in range(21))
DAMPING_VALUES = tuple(10.0**exponent for exponent in (-2.0, -1.0, 0.0, 1.0, 2.0))
RATE_RATIOS = (1.25, 1.75, 3.0, 6.0)
MORPHOLOGY_VALUES = (0.4 + 0.7j, -1.2 + 0.3j, 2.0 - 1.5j)
BASE_FREQUENCIES = (0.35, 1.0, 2.5)
COMMON_GAINS = tuple(
    magnitude * complex(math.cos(phase), math.sin(phase))
    for magnitude, phase in ((0.2, -2.4), (0.75, -0.3), (2.5, 1.1), (8.0, 2.7))
)
ORDER_BOUNDS = (1.0e-6, 1.0)
SOLVER_TOLERANCE = 1.0e-13
CONSISTENCY_TOLERANCE = 1.0e-9
TAMPER_DAMPING_FACTOR = 1.05

FIELDNAMES = (
    "case_id",
    "order_true",
    "damping_true",
    "morphology_real",
    "morphology_imag",
    "frequency_1",
    "frequency_2",
    "frequency_3",
    "frequency_ratio",
    "gain_real",
    "gain_imag",
    "raw_order_abs_error",
    "raw_damping_abs_error",
    "raw_morphology_abs_error",
    "gained_order_abs_error",
    "gained_damping_abs_error",
    "gained_morphology_abs_error",
    "common_gain_order_difference",
    "common_gain_damping_difference",
    "positive_consistent",
    "positive_order_spread",
    "positive_damping_spread",
    "geometric_cross_ratio_abs_error",
    "tamper_consistency_outcome",
    "tamper_detected",
    "tampered_cross_ratio_deviation",
)


def _sha256_json(value: Any) -> str:
    encoded = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _format_scalar(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise RuntimeError("the sweep produced a non-finite value")
        return format(value, ".17g")
    return str(value)


def _complex_record(value: complex) -> dict[str, float]:
    return {"real": value.real, "imag": value.imag}


def _production_source_manifest() -> tuple[str, list[dict[str, str]]]:
    entries = [
        {
            "path": path.relative_to(ROOT).as_posix(),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
        for path in PRODUCTION_SOURCE_PATHS
    ]
    return _sha256_json(entries), entries


def _production_sources_dirty() -> bool:
    completed = subprocess.run(
        [
            "git",
            "status",
            "--porcelain=v1",
            "--untracked-files=all",
            "--",
            *(path.relative_to(ROOT).as_posix() for path in PRODUCTION_SOURCE_PATHS),
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return bool(completed.stdout.strip())


def _current_git_commit() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _is_sha256(value: object) -> bool:
    return (
        type(value) is str
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _is_git_commit(value: object) -> bool:
    return (
        type(value) is str
        and len(value) == 40
        and all(character in "0123456789abcdef" for character in value)
    )


def _read_sha256sums(path: Path) -> dict[str, str]:
    try:
        lines = path.read_text(encoding="ascii").splitlines()
    except (OSError, UnicodeError) as error:
        raise RuntimeError("existing evidence has an unreadable SHA256SUMS") from error
    declared: dict[str, str] = {}
    for line in lines:
        parts = line.split("  ", maxsplit=1)
        if len(parts) != 2:
            raise RuntimeError("existing SHA256SUMS has a malformed entry")
        digest, name = parts
        if not _is_sha256(digest) or name in declared:
            raise RuntimeError("existing SHA256SUMS has an invalid or duplicate entry")
        declared[name] = digest
    return declared


def _recorded_commit_from_existing_evidence(
    output_dir: Path,
    source_manifest_sha256: str,
    source_manifest: list[dict[str, str]],
) -> str | None:
    paths = {name: output_dir / name for name in CANONICAL_ARTIFACT_FILES}
    existence = {name: path.is_file() for name, path in paths.items()}
    if not any(existence.values()):
        return None
    if not all(existence.values()):
        raise RuntimeError("existing evidence is incomplete; archive it before regeneration")

    try:
        stored = json.loads(paths["summary.json"].read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise RuntimeError("existing evidence summary is unreadable") from error
    if type(stored) is not dict or stored.get("schema") != SCHEMA:
        raise RuntimeError("existing evidence summary has an invalid schema")
    unsigned = dict(stored)
    stored_seal = unsigned.pop("summary_sha256", None)
    if not _is_sha256(stored_seal) or stored_seal != _sha256_json(unsigned):
        raise RuntimeError("existing evidence summary has an invalid internal seal")

    declared_artifacts = stored.get("artifact_files")
    if type(declared_artifacts) is not dict or set(declared_artifacts) != {"cases.csv"}:
        raise RuntimeError("existing evidence does not bind its case artifact")
    expected_cases_sha256 = declared_artifacts["cases.csv"]
    if not _is_sha256(expected_cases_sha256):
        raise RuntimeError("existing evidence has an invalid case-artifact hash")
    if hashlib.sha256(paths["cases.csv"].read_bytes()).hexdigest() != expected_cases_sha256:
        raise RuntimeError("existing case artifact hash mismatch")

    declared_sums = _read_sha256sums(paths["SHA256SUMS"])
    if set(declared_sums) != {"cases.csv", "summary.json"}:
        raise RuntimeError("existing SHA256SUMS has missing or unexpected entries")
    for name in ("cases.csv", "summary.json"):
        actual_sha256 = hashlib.sha256(paths[name].read_bytes()).hexdigest()
        if declared_sums[name] != actual_sha256:
            raise RuntimeError("existing SHA256SUMS does not match evidence files")

    provenance = stored.get("provenance")
    if type(provenance) is not dict:
        raise RuntimeError("existing evidence has no provenance record")
    if (
        provenance.get("production_source_manifest_sha256") != source_manifest_sha256
        or provenance.get("production_source_manifest") != source_manifest
    ):
        raise RuntimeError(
            "existing evidence was produced by different production sources; "
            "archive it before regeneration"
        )
    recorded_commit = provenance.get("code_commit")
    if not _is_git_commit(recorded_commit):
        raise RuntimeError("existing evidence has no valid code_commit provenance")
    return recorded_commit


def _frozen_provenance(
    output_dir: Path,
    explicit_git_commit: str | None,
) -> tuple[str, str, list[dict[str, str]]]:
    source_manifest_sha256, source_manifest = _production_source_manifest()
    if explicit_git_commit is not None:
        if not _is_git_commit(explicit_git_commit):
            raise ValueError("explicit git commit must be 40 lowercase hexadecimal characters")
        return explicit_git_commit, source_manifest_sha256, source_manifest
    if _production_sources_dirty():
        raise RuntimeError(
            "refusing frozen evidence generation from dirty production sources; commit code first"
        )
    current_commit = _current_git_commit()
    if not _is_git_commit(current_commit):
        raise RuntimeError("git did not return a valid full commit hash")
    recorded_commit = _recorded_commit_from_existing_evidence(
        output_dir,
        source_manifest_sha256,
        source_manifest,
    )
    return (
        recorded_commit if recorded_commit is not None else current_commit,
        source_manifest_sha256,
        source_manifest,
    )


def _environment_record() -> dict[str, object]:
    return {
        "python_implementation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "platform_machine": platform.machine(),
        "platform_release": platform.release(),
        "platform_system": platform.system(),
        "float_mantissa_bits": sys.float_info.mant_dig,
    }


def _run_case(
    case_id: int,
    order: float,
    damping: float,
    morphology: complex,
    frequency_1: float,
    frequency_ratio: float,
    gain: complex,
) -> dict[str, object]:
    frequency_2 = frequency_1 * frequency_ratio
    frequency_3 = frequency_1 * frequency_ratio**2
    frequencies = (frequency_1, frequency_2, frequency_3)
    responses = tuple(
        harmonic_response(order, damping, morphology, frequency) for frequency in frequencies
    )

    raw = recover_reference_pair(
        responses[:2],
        frequencies[:2],
        order_bounds=ORDER_BOUNDS,
        tolerance=SOLVER_TOLERANCE,
    )
    gained = recover_reference_pair(
        tuple(gain * response for response in responses[:2]),
        frequencies[:2],
        order_bounds=ORDER_BOUNDS,
        tolerance=SOLVER_TOLERANCE,
    )
    positive = evaluate_reference_pair_consistency(
        responses,
        frequencies,
        order_bounds=ORDER_BOUNDS,
        solver_tolerance=SOLVER_TOLERANCE,
        consistency_tolerance=CONSISTENCY_TOLERANCE,
    )
    if not positive.consistent:
        raise RuntimeError(f"exact three-rate case {case_id} failed consistency")

    expected_cross_ratio = 1.0 + frequency_ratio**order
    cross_ratio_error = abs(response_invariant(responses) - expected_cross_ratio)
    tampered_damping = TAMPER_DAMPING_FACTOR * damping
    fractional_frequency_1 = frequency_1**order * complex(
        math.cos(0.5 * math.pi * order),
        math.sin(0.5 * math.pi * order),
    )
    anchored_tampered_morphology = responses[0] * (tampered_damping + fractional_frequency_1)
    tampered_response_3 = harmonic_response(
        order,
        tampered_damping,
        anchored_tampered_morphology,
        frequency_3,
    )
    tampered_responses = (responses[0], responses[1], tampered_response_3)
    tampered_cross_ratio_deviation = abs(
        response_invariant(tampered_responses) - expected_cross_ratio
    )
    try:
        tampered = evaluate_reference_pair_consistency(
            tampered_responses,
            frequencies,
            order_bounds=ORDER_BOUNDS,
            solver_tolerance=SOLVER_TOLERANCE,
            consistency_tolerance=CONSISTENCY_TOLERANCE,
        )
    except ValueError:
        tamper_outcome = "REJECTED_PHYSICAL_PAIR"
        tamper_detected = True
    else:
        tamper_outcome = "INCONSISTENT" if not tampered.consistent else "MISSED"
        tamper_detected = not tampered.consistent
    if not tamper_detected or tampered_cross_ratio_deviation <= 1.0e-8:
        raise RuntimeError(f"deterministic third-rate tamper was not detected in case {case_id}")

    expected_gained_morphology = gain * morphology
    return {
        "case_id": case_id,
        "order_true": order,
        "damping_true": damping,
        "morphology_real": morphology.real,
        "morphology_imag": morphology.imag,
        "frequency_1": frequency_1,
        "frequency_2": frequency_2,
        "frequency_3": frequency_3,
        "frequency_ratio": frequency_ratio,
        "gain_real": gain.real,
        "gain_imag": gain.imag,
        "raw_order_abs_error": abs(raw.order - order),
        "raw_damping_abs_error": abs(raw.damping - damping),
        "raw_morphology_abs_error": abs(raw.morphology - morphology),
        "gained_order_abs_error": abs(gained.order - order),
        "gained_damping_abs_error": abs(gained.damping - damping),
        "gained_morphology_abs_error": abs(gained.morphology - expected_gained_morphology),
        "common_gain_order_difference": abs(gained.order - raw.order),
        "common_gain_damping_difference": abs(gained.damping - raw.damping),
        "positive_consistent": positive.consistent,
        "positive_order_spread": positive.maximum_order_spread,
        "positive_damping_spread": positive.maximum_damping_spread,
        "geometric_cross_ratio_abs_error": cross_ratio_error,
        "tamper_consistency_outcome": tamper_outcome,
        "tamper_detected": tamper_detected,
        "tampered_cross_ratio_deviation": tampered_cross_ratio_deviation,
    }


def _case_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    case_id = 0
    for order in ORDER_VALUES:
        for damping in DAMPING_VALUES:
            for ratio in RATE_RATIOS:
                for morphology_index, morphology in enumerate(MORPHOLOGY_VALUES):
                    rows.append(
                        _run_case(
                            case_id,
                            order,
                            damping,
                            morphology,
                            BASE_FREQUENCIES[morphology_index],
                            ratio,
                            COMMON_GAINS[case_id % len(COMMON_GAINS)],
                        )
                    )
                    case_id += 1
    return rows


def _maximum(rows: list[dict[str, object]], field: str) -> float:
    return max(float(row[field]) for row in rows)


def _render_cases(rows: list[dict[str, object]]) -> str:
    handle = io.StringIO(newline="")
    writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({name: _format_scalar(row[name]) for name in FIELDNAMES})
    return handle.getvalue()


def run_two_rate_recovery_sweep(
    output_dir: Path = DEFAULT_OUTPUT,
    *,
    git_commit: str | None = None,
) -> dict[str, Any]:
    """Execute the grid, persist replayable evidence, and return its summary."""

    bound_git_commit, source_manifest_sha256, source_manifest = _frozen_provenance(
        output_dir,
        git_commit,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = _case_rows()
    expected_cases = (
        len(ORDER_VALUES) * len(DAMPING_VALUES) * len(RATE_RATIOS) * len(MORPHOLOGY_VALUES)
    )
    if len(rows) != expected_cases or expected_cases < 1200:
        raise RuntimeError("the deterministic grid did not meet the declared case count")

    cases_text = _render_cases(rows)
    repeated_cases_text = _render_cases(_case_rows())
    case_serialization_byte_identical = cases_text == repeated_cases_text
    if not case_serialization_byte_identical:
        raise RuntimeError("the repeated deterministic grid was not byte-identical")

    if _production_source_manifest()[0] != source_manifest_sha256:
        raise RuntimeError("production sources changed during the sweep")
    cases_path = output_dir / "cases.csv"
    cases_path.write_text(cases_text, encoding="utf-8", newline="\n")
    cases_sha256 = hashlib.sha256(cases_path.read_bytes()).hexdigest()
    tamper_counts: dict[str, int] = {}
    for row in rows:
        outcome = str(row["tamper_consistency_outcome"])
        tamper_counts[outcome] = tamper_counts.get(outcome, 0) + 1

    summary: dict[str, Any] = {
        "schema": SCHEMA,
        "claim_boundary": {
            "analytical_proof_replaced": False,
            "noisy_or_real_ecg_validated": False,
            "purpose": "floating-point audit of exact-model global inverse identities",
        },
        "design": {
            "generator": (
                "fixed Cartesian product of order x damping x rate ratio x "
                "paired base-frequency/morphology condition; no PRNG"
            ),
            "seed": None,
            "case_count": len(rows),
            "order_values": list(ORDER_VALUES),
            "damping_values": list(DAMPING_VALUES),
            "frequency_ratios": list(RATE_RATIOS),
            "base_frequency_morphology_pairs": [
                {
                    "base_frequency": frequency,
                    "morphology": _complex_record(morphology),
                }
                for frequency, morphology in zip(BASE_FREQUENCIES, MORPHOLOGY_VALUES, strict=True)
            ],
            "pairing_rule": "paired by the shared condition index, not a Cartesian cross",
            "common_gains": [_complex_record(value) for value in COMMON_GAINS],
            "order_bounds": list(ORDER_BOUNDS),
            "solver_tolerance": SOLVER_TOLERANCE,
            "consistency_tolerance": CONSISTENCY_TOLERANCE,
            "tamper_third_pair_damping_factor": TAMPER_DAMPING_FACTOR,
            "tamper_anchor": "first response held exact; only third response replaced",
        },
        "determinism": {
            "case_serialization_repeats": 2,
            "case_serialization_byte_identical": case_serialization_byte_identical,
            "fixed_grid_no_prng": True,
        },
        "provenance": {
            "code_commit": bound_git_commit,
            "production_source_manifest_sha256": source_manifest_sha256,
            "production_source_manifest": source_manifest,
            "environment": _environment_record(),
        },
        "outcomes": {
            "two_rate_global_inverse": {
                "successful_cases": len(rows),
                "failed_cases": 0,
                "maximum_order_absolute_error": _maximum(rows, "raw_order_abs_error"),
                "maximum_damping_absolute_error": _maximum(rows, "raw_damping_abs_error"),
                "maximum_morphology_absolute_error": _maximum(rows, "raw_morphology_abs_error"),
            },
            "common_complex_gain_invariance": {
                "successful_cases": len(rows),
                "failed_cases": 0,
                "maximum_order_difference": _maximum(rows, "common_gain_order_difference"),
                "maximum_damping_difference": _maximum(rows, "common_gain_damping_difference"),
                "maximum_gained_morphology_absolute_error": _maximum(
                    rows, "gained_morphology_abs_error"
                ),
            },
            "three_rate_positive_consistency": {
                "consistent_cases": sum(bool(row["positive_consistent"]) for row in rows),
                "inconsistent_cases": sum(not bool(row["positive_consistent"]) for row in rows),
                "maximum_order_spread": _maximum(rows, "positive_order_spread"),
                "maximum_damping_spread": _maximum(rows, "positive_damping_spread"),
                "maximum_geometric_cross_ratio_absolute_error": _maximum(
                    rows, "geometric_cross_ratio_abs_error"
                ),
            },
            "deterministic_third_rate_tamper": {
                "detected_cases": sum(bool(row["tamper_detected"]) for row in rows),
                "missed_cases": sum(not bool(row["tamper_detected"]) for row in rows),
                "consistency_outcome_counts": dict(sorted(tamper_counts.items())),
                "minimum_tampered_cross_ratio_deviation": min(
                    float(row["tampered_cross_ratio_deviation"]) for row in rows
                ),
            },
        },
        "artifact_files": {"cases.csv": cases_sha256},
    }
    summary["summary_sha256"] = _sha256_json(summary)
    summary_path = output_dir / "summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, ensure_ascii=True, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    summary_file_sha256 = hashlib.sha256(summary_path.read_bytes()).hexdigest()
    (output_dir / "SHA256SUMS").write_text(
        f"{cases_sha256}  cases.csv\n{summary_file_sha256}  summary.json\n",
        encoding="ascii",
        newline="\n",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    summary = run_two_rate_recovery_sweep(args.output_dir.resolve())
    print(json.dumps(summary["outcomes"], indent=2, sort_keys=True))
    print(f"summary_sha256={summary['summary_sha256']}")


if __name__ == "__main__":
    main()
