import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "literature" / "manifest.csv"
LOCK = ROOT / "literature" / "download-lock.json"


def load_manifest() -> list[dict[str, str]]:
    with MANIFEST.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def test_literature_manifest_has_unique_reproducible_entries() -> None:
    rows = load_manifest()
    identifiers = [row["id"] for row in rows]
    filenames = [row["local_filename"] for row in rows if row["local_filename"]]

    assert len(rows) >= 20
    assert len(identifiers) == len(set(identifiers))
    assert len(filenames) == len(set(filenames))
    assert all(row["landing_url"].startswith("https://") for row in rows)

    for row in rows:
        if row["download"] == "1":
            assert row["pdf_url"].startswith("https://")
            assert row["local_filename"].endswith(".pdf")


def test_literature_lock_matches_selected_manifest() -> None:
    rows = load_manifest()
    selected = {row["local_filename"] for row in rows if row["download"] == "1"}
    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    files = lock["files"]

    assert {item["file"] for item in files} == selected
    assert lock["total_bytes"] == sum(item["bytes"] for item in files)
    assert all(0 < item["bytes"] <= lock["max_file_bytes"] for item in files)
    assert all(item["pages"] > 0 for item in files)
    assert all(re.fullmatch(r"[0-9a-f]{64}", item["sha256"]) for item in files)
