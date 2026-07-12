from __future__ import annotations

import hashlib
from pathlib import Path

from scripts.replay_fantasia_wls import BUNDLE_FILES, validate_bundle_checksums


def _write_sealed_bundle(root: Path) -> None:
    lines = []
    for index, name in enumerate(BUNDLE_FILES):
        payload = f"fixture-{index}-{name}\n".encode()
        (root / name).write_bytes(payload)
        lines.append(f"{hashlib.sha256(payload).hexdigest()}  {name}")
    (root / "SHA256SUMS").write_text("\n".join(lines) + "\n", encoding="ascii")


def test_bundle_checksum_allowlist_accepts_exact_fixture(tmp_path: Path) -> None:
    _write_sealed_bundle(tmp_path)

    assert validate_bundle_checksums(tmp_path) == []


def test_bundle_checksum_allowlist_rejects_tamper_and_extra_entry(tmp_path: Path) -> None:
    _write_sealed_bundle(tmp_path)
    (tmp_path / BUNDLE_FILES[0]).write_bytes(b"tampered")
    with (tmp_path / "SHA256SUMS").open("a", encoding="ascii") as stream:
        stream.write(f"{'0' * 64}  unexpected.bin\n")

    reasons = {row["reason"] for row in validate_bundle_checksums(tmp_path)}

    assert "bundle_file_sha256_mismatch" in reasons
    assert "SHA256SUMS_allowlist_mismatch" in reasons


def test_bundle_checksum_allowlist_rejects_duplicate_entry(tmp_path: Path) -> None:
    _write_sealed_bundle(tmp_path)
    first = (tmp_path / "SHA256SUMS").read_text(encoding="ascii").splitlines()[0]
    with (tmp_path / "SHA256SUMS").open("a", encoding="ascii") as stream:
        stream.write(first + "\n")

    reasons = {row["reason"] for row in validate_bundle_checksums(tmp_path)}

    assert "duplicate_SHA256SUMS_entry" in reasons
