"""Replay the serialized Fantasia WLS execution records from external data."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import wfdb

from fo_ekf.fantasia_pilot import (
    FANTASIA_V1_SHA256,
    FantasiaPilotConfig,
    canonical_json_hash,
    fantasia_implementation_manifest,
    fantasia_runtime_environment,
    sha256_file,
)
from fo_ekf.wls_execution import (
    bounded_sample_indices,
    replay_wls_execution_record,
    validate_execution_manifest_binding,
)

BUNDLE_FILES = (
    "design_lock.json",
    "window_manifest.csv.gz",
    "wls_coefficients.csv.gz",
    "wls_execution_records.jsonl.gz",
    "heldout_metrics.csv",
    "model_fits.json",
    "source_inventory.json",
    "r2_status.json",
    "summary.json",
)


def validate_bundle_checksums(bundle_dir: Path) -> list[dict[str, object]]:
    """Validate the exact canonical bundle allowlist before parsing evidence."""

    sums_path = bundle_dir / "SHA256SUMS"
    if not sums_path.is_file():
        return [{"reason": "missing_SHA256SUMS"}]
    declared: dict[str, str] = {}
    failures: list[dict[str, object]] = []
    for line in sums_path.read_text(encoding="ascii").splitlines():
        try:
            digest, name = line.split("  ", maxsplit=1)
        except ValueError:
            failures.append({"reason": "malformed_SHA256SUMS_line"})
            continue
        if name in declared:
            failures.append({"reason": "duplicate_SHA256SUMS_entry", "file": name})
        declared[name] = digest
    if set(declared) != set(BUNDLE_FILES):
        failures.append(
            {
                "reason": "SHA256SUMS_allowlist_mismatch",
                "missing": sorted(set(BUNDLE_FILES) - set(declared)),
                "unexpected": sorted(set(declared) - set(BUNDLE_FILES)),
            }
        )
    for name in BUNDLE_FILES:
        path = bundle_dir / name
        if not path.is_file():
            failures.append({"reason": "bundle_file_missing", "file": name})
        elif name in declared and hashlib.sha256(path.read_bytes()).hexdigest() != declared[name]:
            failures.append({"reason": "bundle_file_sha256_mismatch", "file": name})
    return failures


def _default_data_dir() -> Path:
    raw = os.environ.get("ECG_DATA_DIR")
    return Path(raw) if raw else Path(__file__).resolve().parents[2] / "PUBLIC" / "fantasia"


def _default_bundle_dir() -> Path:
    raw = os.environ.get("ECG_OUTPUT_DIR")
    if raw:
        return Path(raw) / "fantasia_pilot_candidate"
    return (
        Path(__file__).resolve().parents[1] / "research" / "results" / "fantasia_pilot_2026-07-12"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=_default_data_dir())
    parser.add_argument("--bundle-dir", type=Path, default=_default_bundle_dir())
    arguments = parser.parse_args()
    data_dir = arguments.data_dir.resolve()
    bundle_dir = arguments.bundle_dir.resolve()

    checksum_failures = validate_bundle_checksums(bundle_dir)
    if checksum_failures:
        print(
            json.dumps(
                {
                    "status": "FAIL",
                    "failure_count": len(checksum_failures),
                    "failures": checksum_failures,
                    "certification_status": "NOT_CERTIFIABLE",
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 1

    source_inventory = json.loads((bundle_dir / "source_inventory.json").read_text("utf-8"))
    stored_record_hashes = source_inventory["record_sha256"]
    design_lock = json.loads((bundle_dir / "design_lock.json").read_text("utf-8"))
    summary = json.loads((bundle_dir / "summary.json").read_text("utf-8"))
    protocol_sha256 = design_lock["config_sha256"]
    with gzip.open(
        bundle_dir / "window_manifest.csv.gz", "rt", encoding="utf-8", newline=""
    ) as stream:
        manifest_rows = list(csv.DictReader(stream))
    with gzip.open(bundle_dir / "wls_execution_records.jsonl.gz", "rt", encoding="utf-8") as stream:
        execution_rows = [json.loads(line) for line in stream if line.strip()]
    config = FantasiaPilotConfig()
    if {row["record"] for row in execution_rows} != set(config.records):
        raise ValueError("serialized execution records do not match the locked record set")

    failures: list[dict[str, object]] = []
    current_config_sha256 = config.lock_hash()
    current_implementation = fantasia_implementation_manifest()
    current_environment = fantasia_runtime_environment()
    if protocol_sha256 != current_config_sha256 or summary.get("config_sha256") != protocol_sha256:
        failures.append({"reason": "config_lock_mismatch"})
    if (
        design_lock.get("implementation_manifest") != current_implementation
        or summary.get("implementation_sha256") != current_implementation["implementation_sha256"]
    ):
        failures.append({"reason": "production_implementation_manifest_mismatch"})
    if (
        design_lock.get("runtime_environment") != current_environment
        or summary.get("runtime_environment") != current_environment
    ):
        failures.append({"reason": "runtime_environment_mismatch"})
    dataset = source_inventory.get("dataset")
    if not isinstance(dataset, dict) or (
        dataset.get("name"),
        dataset.get("version"),
        dataset.get("doi"),
        dataset.get("file_license"),
        dataset.get("selected_files_match_upstream_sha256_manifest"),
    ) != (
        "Fantasia Database",
        "1.0.0",
        "10.13026/C2RG61",
        "Open Data Commons Attribution License v1.0",
        True,
    ):
        failures.append({"reason": "dataset_provenance_mismatch"})
    inventory_hashes = {
        f"{row['record']}.{row['extension']}": row["sha256"]
        for row in source_inventory.get("files", [])
        if isinstance(row, dict)
        and isinstance(row.get("record"), str)
        and isinstance(row.get("extension"), str)
        and isinstance(row.get("sha256"), str)
    }
    if inventory_hashes != FANTASIA_V1_SHA256:
        failures.append({"reason": "source_inventory_not_official_fantasia_v1"})
    if failures:
        print(
            json.dumps(
                {
                    "status": "FAIL",
                    "execution_record_count": len(execution_rows),
                    "manifest_record_count": len(manifest_rows),
                    "matched_count": 0,
                    "failure_count": len(failures),
                    "failures": failures,
                    "certification_status": "NOT_CERTIFIABLE",
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 1
    manifest_by_key: dict[tuple[str, str, int], dict[str, str]] = {}
    for row in manifest_rows:
        key = (row["record"], row["split"], int(row["window_index"]))
        if key in manifest_by_key:
            failures.append({"key": key, "reason": "duplicate_manifest_key"})
        manifest_by_key[key] = row
    execution_by_key: dict[tuple[str, str, int], dict[str, object]] = {}
    for row in execution_rows:
        key = (row["record"], row["split"], int(row["window_index"]))
        if key in execution_by_key:
            failures.append({"key": key, "reason": "duplicate_execution_key"})
        execution_by_key[key] = row
    missing_manifest = sorted(set(execution_by_key) - set(manifest_by_key))
    missing_execution = sorted(set(manifest_by_key) - set(execution_by_key))
    if missing_manifest:
        failures.append({"reason": "execution_without_manifest", "keys": missing_manifest})
    if missing_execution:
        failures.append({"reason": "manifest_without_execution", "keys": missing_execution})

    matched = 0
    for record_name in config.records:
        file_hashes = {
            extension: sha256_file(data_dir / f"{record_name}.{extension}")
            for extension in ("hea", "dat", "ecg")
        }
        if any(
            file_hashes[extension] != FANTASIA_V1_SHA256[f"{record_name}.{extension}"]
            for extension in ("hea", "dat", "ecg")
        ):
            failures.append({"record": record_name, "reason": "upstream_file_sha256_mismatch"})
            continue
        record_sha256 = canonical_json_hash({"record": record_name, "files": file_hashes})
        if record_sha256 != stored_record_hashes[record_name]:
            failures.append({"record": record_name, "reason": "external_record_hash_mismatch"})
            continue

        header = wfdb.rdheader(str(data_dir / record_name))
        channel_index = tuple(header.sig_name).index("ECG")
        record = wfdb.rdrecord(
            str(data_dir / record_name), channels=[channel_index], physical=False
        )
        digital = np.asarray(record.d_signal[:, 0], dtype=np.int64)
        for row in (item for item in execution_rows if item["record"] == record_name):
            key = (row["record"], row["split"], int(row["window_index"]))
            manifest_row = manifest_by_key.get(key)
            if manifest_row is None:
                continue
            binding_failures = validate_execution_manifest_binding(
                row,
                manifest_row,
                external_record_sha256=record_sha256,
                protocol_sha256=protocol_sha256,
            )
            if binding_failures:
                failures.append(
                    {
                        "record": record_name,
                        "split": row["split"],
                        "window_index": row["window_index"],
                        "reason": "cross_binding_mismatch",
                        "fields": binding_failures,
                    }
                )
                continue
            document = row["execution"]
            window = document["window"]
            selected_indices = bounded_sample_indices(
                window["start_sample_inclusive"],
                window["end_sample_exclusive"],
            )
            payload = digital[selected_indices]
            replay = replay_wls_execution_record(document, payload)
            if replay["status"] == "MATCH":
                matched += 1
            else:
                failures.append(
                    {
                        "record": record_name,
                        "split": row["split"],
                        "window_index": row["window_index"],
                        "replay": replay,
                    }
                )

    result = {
        "status": "PASS" if not failures else "FAIL",
        "execution_record_count": len(execution_rows),
        "manifest_record_count": len(manifest_rows),
        "matched_count": matched,
        "failure_count": len(failures),
        "failures": failures,
        "certification_status": "NOT_CERTIFIABLE",
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
