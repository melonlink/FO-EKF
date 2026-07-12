"""Run the deterministic exact-sample WLS -> R2 -> T1 synthetic closure.

The production synthetic-case builder supplies a directly constructed
steady-state harmonic signal and declared generator bounds.  Every scientific
relation is recomputed by production certificate code, and the emitted summary
records both positive replay and tamper controls.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
for path in (str(ROOT), str(SRC)):
    if path not in sys.path:
        sys.path.insert(0, path)

from fo_ekf.r2_bound_certificate import replay_r2_bound_certificate  # noqa: E402
from fo_ekf.r2_protocol import validate_r2_protocol  # noqa: E402
from fo_ekf.r2_t1_bundle import (  # noqa: E402
    R2T1BundleRelation,
    certify_r2_t1_bundle,
    replay_r2_t1_bundle,
)
from fo_ekf.synthetic_certificate_case import build_synthetic_certificate_case  # noqa: E402
from fo_ekf.wls_execution import replay_wls_execution_record  # noqa: E402
from fo_ekf.wls_interval_certificate import replay_wls_interval_certificate  # noqa: E402

DEFAULT_OUTPUT = ROOT / "research" / "results" / "synthetic_certificate_chain_2026-07-12"


def _sha256_json(value: Any) -> str:
    encoded = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _production_source_manifest() -> tuple[str, list[dict[str, str]]]:
    paths = [
        ROOT / "pyproject.toml",
        ROOT / "requirements-repro-lock.txt",
        ROOT / "config" / "r2_protocol_template.toml",
        ROOT / "scripts" / "run_synthetic_certificate_chain.py",
        *sorted((ROOT / "src" / "fo_ekf").glob("*.py")),
    ]
    entries = [
        {
            "path": path.relative_to(ROOT).as_posix(),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
        for path in paths
    ]
    return _sha256_json(entries), entries


def _frozen_provenance(
    output_dir: Path,
    explicit_git_commit: str | None,
) -> tuple[str, str, list[dict[str, str]]]:
    source_sha256, source_entries = _production_source_manifest()
    if explicit_git_commit is not None:
        return explicit_git_commit, source_sha256, source_entries
    status = subprocess.run(
        [
            "git",
            "status",
            "--porcelain=v1",
            "--untracked-files=all",
            "--",
            "pyproject.toml",
            "requirements-repro-lock.txt",
            "config/r2_protocol_template.toml",
            "scripts/run_synthetic_certificate_chain.py",
            "src/fo_ekf",
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    if status.stdout.strip():
        raise RuntimeError(
            "refusing frozen evidence generation from dirty production sources; commit code first"
        )
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    current_commit = completed.stdout.strip()
    existing_summary = output_dir / "summary.json"
    if existing_summary.is_file():
        stored = json.loads(existing_summary.read_text(encoding="utf-8"))
        stored_unsigned = dict(stored)
        stored_summary_sha256 = stored_unsigned.pop("summary_sha256", None)
        if stored_summary_sha256 != _sha256_json(stored_unsigned):
            raise RuntimeError("existing evidence summary has an invalid internal seal")
        declared_artifacts = stored.get("artifact_files")
        if type(declared_artifacts) is not dict or set(declared_artifacts) != {"artifacts.json"}:
            raise RuntimeError("existing evidence does not bind its replay artifacts")
        for name, expected_sha256 in declared_artifacts.items():
            artifact_path = output_dir / name
            if (
                type(name) is not str
                or type(expected_sha256) is not str
                or not artifact_path.is_file()
                or hashlib.sha256(artifact_path.read_bytes()).hexdigest() != expected_sha256
            ):
                raise RuntimeError("existing replay artifact hash mismatch")
        sums_path = output_dir / "SHA256SUMS"
        if not sums_path.is_file():
            raise RuntimeError("existing evidence has no SHA256SUMS")
        declared_sums = {}
        for line in sums_path.read_text(encoding="ascii").splitlines():
            digest, name = line.split("  ", maxsplit=1)
            if name in declared_sums:
                raise RuntimeError("existing SHA256SUMS contains a duplicate entry")
            declared_sums[name] = digest
        if set(declared_sums) != {"artifacts.json", "summary.json"}:
            raise RuntimeError("existing SHA256SUMS has missing or unexpected entries")
        for name in ("artifacts.json", "summary.json"):
            path = output_dir / name
            if (
                name not in declared_sums
                or not path.is_file()
                or hashlib.sha256(path.read_bytes()).hexdigest() != declared_sums[name]
            ):
                raise RuntimeError("existing SHA256SUMS does not match evidence files")
        stored_scope = stored.get("input_scope", {})
        if stored_scope.get("production_source_manifest_sha256") != source_sha256:
            raise RuntimeError(
                "existing evidence was produced by different production sources; "
                "archive it before regeneration"
            )
        recorded_commit = stored_scope.get("code_commit")
        if not isinstance(recorded_commit, str):
            raise RuntimeError("existing evidence has no code_commit provenance")
        return recorded_commit, source_sha256, source_entries
    return current_commit, source_sha256, source_entries


def run_synthetic_chain(
    output_dir: Path,
    *,
    git_commit: str | None = None,
) -> dict[str, Any]:
    """Execute, verify, adversarially perturb, and persist one small chain."""

    bound_git_commit, source_manifest_sha256, source_manifest = _frozen_provenance(
        output_dir,
        git_commit,
    )
    case = build_synthetic_certificate_case(git_commit=bound_git_commit)
    problem = case.problem
    box = case.box
    baseline_json = case.baseline_certificate_json
    protocol = case.protocol
    links = case.links
    link = links[0]
    execution = link.wls_execution_record
    interval = link.wls_interval_certificate
    if not isinstance(execution, dict) or not isinstance(interval, dict):
        raise RuntimeError("synthetic material did not provide complete WLS evidence")
    digital = np.asarray(link.digital_samples, dtype=np.int64)

    execution_replay = replay_wls_execution_record(execution, digital)
    interval_replay = replay_wls_interval_certificate(interval, execution, digital)
    component_replay = replay_r2_bound_certificate(link.request, link.json_text())
    protocol_replay = validate_r2_protocol(protocol)
    bundle = certify_r2_t1_bundle(problem, box, protocol, baseline_json, links)
    bundle_replay = replay_r2_t1_bundle(
        problem,
        box,
        protocol,
        baseline_json,
        links,
        bundle,
    )
    if execution_replay.get("status") != "MATCH":
        raise RuntimeError(f"WLS execution did not replay: {execution_replay}")
    if (
        interval_replay.get("status") != "MATCH"
        or interval_replay.get("relation") != "CERTIFIED_INTERVAL"
    ):
        raise RuntimeError(f"exact interval did not replay: {interval_replay}")
    if not component_replay.valid or component_replay.relation.value != "CERTIFIED_BOUND":
        raise RuntimeError(f"R2 component did not replay: {component_replay}")
    if protocol_replay.status.value != "PASS_DETERMINISTIC":
        raise RuntimeError(f"protocol did not pass: {protocol_replay.reasons}")
    if bundle.relation is not R2T1BundleRelation.PRESERVED_ROBUST or not bundle_replay.valid:
        raise RuntimeError(f"bundle did not close: {bundle.reason}/{bundle_replay.reason}")

    component_document = json.loads(link.json_text())
    component_sampling = Fraction(
        component_document["proof"]["sampling"]["coefficient_radius_upper"]
    )
    exact_row = interval["response_disks"][0]
    exact_dyadic = exact_row["absolute_error_upper"]["dyadic"]
    exact_floor = Fraction(int(exact_dyadic["numerator"]), int(exact_dyadic["denominator"]))
    window_row = protocol["window_result"][0]
    declared_sampling = Fraction.from_float(window_row["radius_sampling"])
    if declared_sampling < component_sampling + exact_floor:
        raise RuntimeError("radius_sampling is below component plus exact-functional floor")
    known_response = problem.harmonics[0].measured_responses[0]
    published_center = bundle.response_disks[0].center
    published_radius = bundle.response_disks[0].radius
    center_absolute_error = abs(published_center - known_response)
    disk_contains_known_response = center_absolute_error <= published_radius
    radius_to_response_ratio = published_radius / abs(known_response)
    if not disk_contains_known_response:
        raise RuntimeError("published response disk does not contain the known synthetic response")
    if not radius_to_response_ratio < 0.01:
        raise RuntimeError("synthetic response disk is too wide to be an informative positive case")

    changed_samples = list(link.digital_samples or ())
    changed_samples[0] += 1
    payload_tamper = replace(link, digital_samples=tuple(changed_samples))
    payload_tamper_result = certify_r2_t1_bundle(
        problem,
        box,
        protocol,
        baseline_json,
        (payload_tamper,),
    )
    changed_interval = deepcopy(interval)
    changed_interval["response_disks"][0]["absolute_error_upper"]["float_upper"] *= 2.0
    interval_tamper = replace(link, wls_interval_certificate=changed_interval)
    interval_tamper_result = certify_r2_t1_bundle(
        problem,
        box,
        protocol,
        baseline_json,
        (interval_tamper,),
    )
    if any(
        result.relation is not R2T1BundleRelation.NOT_CERTIFIABLE
        for result in (payload_tamper_result, interval_tamper_result)
    ):
        raise RuntimeError("a tamper control did not fail closed")

    bundle_document = json.loads(bundle.to_json())
    artifact_document = {
        "schema": "fo-ekf.synthetic-exact-certificate-artifacts.v1",
        "code_commit": bound_git_commit,
        "production_source_manifest_sha256": source_manifest_sha256,
        "production_source_manifest": source_manifest,
        "digital_samples_int64": [int(value) for value in digital],
        "protocol": protocol,
        "baseline_t1_certificate": json.loads(baseline_json),
        "wls_execution": execution,
        "wls_interval_certificate": interval,
        "r2_component_certificate": component_document,
        "r2_t1_bundle_certificate": bundle_document,
    }
    artifact_text = (
        json.dumps(
            artifact_document,
            indent=2,
            sort_keys=True,
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n"
    )
    artifact_sha256 = hashlib.sha256(artifact_text.encode("utf-8")).hexdigest()
    summary: dict[str, Any] = {
        "schema": "fo-ekf.synthetic-exact-certificate-chain.v1",
        "run_id": "synthetic-exact-sample-chain-2026-07-12",
        "input_scope": {
            "kind": "deterministic_synthetic_integer_payload",
            "raw_ecg_used": False,
            "downloads_performed": False,
            "code_commit": bound_git_commit,
            "production_source_manifest_sha256": source_manifest_sha256,
            "sample_count": int(digital.size),
            "sample_times_s": list(link.request.sample_times_s),
            "sample_indices": list(link.request.sample_indices or ()),
            "window_start_sample": link.request.window_start_sample,
            "window_end_sample_exclusive": link.request.window_end_sample_exclusive,
            "phase_operator": "2*pi*N*(sample_index-start)/(end-start)",
            "generator_contract": {
                "kind": "direct_steady_state_retained_harmonic",
                "history_bound": 0.0,
                "dwell_bound": 0.0,
                "timestamp_error_s": 0.0,
                "interpolation_error": 0.0,
                "anti_alias_error": 0.0,
                "measurement_noise": 0.0,
                "model_residual": 0.0,
                "nonzero_terms": [
                    "ADC quantization",
                    "exact-functional numerical enclosure",
                ],
                "scope": "synthetic construction only; zeros are not estimated from ECG",
            },
            "digital_sample_sha256": execution["payload_hashes"]["digital_sample_sha256"],
        },
        "stages": {
            "wls_execution": execution_replay,
            "exact_interval": interval_replay,
            "r2_component": {
                "valid": component_replay.valid,
                "relation": component_replay.relation.value,
                "certificate_sha256": component_document["certificate_sha256"],
            },
            "protocol": {
                "status": protocol_replay.status.value,
                "canonical_sha256": protocol_replay.canonical_sha256,
            },
            "r2_t1_bundle": {
                "relation": bundle.relation.value,
                "reason": bundle.reason,
                "certificate_sha256": bundle.certificate_sha256,
                "replay_valid": bundle_replay.valid,
                "replay_reason": bundle_replay.reason,
            },
        },
        "numeric_radius_contract": {
            "component_sampling_radius_floor": str(component_sampling),
            "exact_functional_radius_floor": str(exact_floor),
            "declared_radius_sampling": str(declared_sampling),
            "composition": "component_sampling_plus_full_exact_functional_floor",
            "solve_only_floor_added_again": False,
        },
        "coverage": {
            "mathematical_exact_sample_wls_coverage": bundle_document["proof"][
                "mathematical_exact_sample_wls_coverage"
            ],
            "uncovered_numeric_term": bundle_document["proof"]["uncovered_numeric_term"],
            "interval_claim_scope": interval["claim_scope"],
            "non_numeric_interval_exclusions": interval["exclusions"],
        },
        "response_disk": {
            "known_model_response_real": known_response.real,
            "known_model_response_imag": known_response.imag,
            "center_real": published_center.real,
            "center_imag": published_center.imag,
            "center_absolute_error": center_absolute_error,
            "radius": published_radius,
            "disk_contains_known_response": disk_contains_known_response,
            "radius_to_response_magnitude_ratio": radius_to_response_ratio,
            "known_parameters": {
                "alpha": 0.7,
                "lambda": 0.5,
                "q_real": 1.0,
                "q_imag": -0.2,
                "angular_frequency_rad_s": 2.0 * np.pi,
            },
        },
        "tamper_controls": {
            "digital_payload_change": {
                "relation": payload_tamper_result.relation.value,
                "reason": payload_tamper_result.reason,
            },
            "interval_radius_change_without_valid_seal": {
                "relation": interval_tamper_result.relation.value,
                "reason": interval_tamper_result.reason,
            },
        },
        "claim_boundary": (
            "This closes the frozen synthetic integer-sample functional and declared R2/T1 "
            "bounds only; it is not real-ECG validation or a cardiac digital-twin claim."
        ),
        "artifact_files": {
            "artifacts.json": artifact_sha256,
        },
    }
    summary["summary_sha256"] = _sha256_json(summary)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_text = (
        json.dumps(summary, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n"
    )
    summary_sha256 = hashlib.sha256(summary_text.encode("utf-8")).hexdigest()
    (output_dir / "artifacts.json").write_bytes(artifact_text.encode("utf-8"))
    (output_dir / "summary.json").write_bytes(summary_text.encode("utf-8"))
    (output_dir / "SHA256SUMS").write_bytes(
        f"{artifact_sha256}  artifacts.json\n{summary_sha256}  summary.json\n".encode("ascii")
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    summary = run_synthetic_chain(arguments.output.resolve())
    print(json.dumps(summary["stages"], indent=2, sort_keys=True))
    print(f"summary: {arguments.output.resolve() / 'summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
