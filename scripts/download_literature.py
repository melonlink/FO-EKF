"""Download the curated open-access literature set with a hard 100 MB limit."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import ssl
import sys
import urllib.request
from pathlib import Path
from typing import BinaryIO

import certifi

# User approval is required above 100 decimal MB, not 100 MiB.
MAX_FILE_BYTES = 100_000_000
CHUNK_BYTES = 1024 * 1024
USER_AGENT = "FO-EKF-literature-audit/0.1 (mailto:melonlink@tsnu.edu.cn)"
SSL_CONTEXT = ssl.create_default_context(cafile=certifi.where())


def sha256_stream(stream: BinaryIO) -> tuple[str, int]:
    digest = hashlib.sha256()
    total = 0
    while chunk := stream.read(CHUNK_BYTES):
        digest.update(chunk)
        total += len(chunk)
    return digest.hexdigest(), total


def sha256_file(path: Path) -> tuple[str, int]:
    with path.open("rb") as stream:
        return sha256_stream(stream)


def is_pdf(path: Path) -> bool:
    with path.open("rb") as stream:
        if stream.read(5) != b"%PDF-":
            return False
        tail_size = min(path.stat().st_size, 64 * 1024)
        stream.seek(-tail_size, 2)
        return b"%%EOF" in stream.read()


def download_pdf(url: str, destination: Path) -> tuple[str, int]:
    temporary = destination.with_suffix(destination.suffix + ".part")
    temporary.unlink(missing_ok=True)
    digest = hashlib.sha256()
    total = 0
    expected_bytes: int | None = None

    try:
        with temporary.open("wb") as stream:
            for _attempt in range(8):
                headers = {"Accept-Encoding": "identity", "User-Agent": USER_AGENT}
                if total:
                    headers["Range"] = f"bytes={total}-"
                request = urllib.request.Request(url, headers=headers)

                with urllib.request.urlopen(  # noqa: S310
                    request, timeout=90, context=SSL_CONTEXT
                ) as response:
                    status = response.status
                    content_length = response.headers.get("Content-Length")
                    content_range = response.headers.get("Content-Range")

                    if total and status != 206:
                        raise RuntimeError("server does not support resuming a truncated response")
                    if content_range:
                        unit, range_value = content_range.split(" ", 1)
                        byte_range, announced_total = range_value.split("/", 1)
                        range_start = int(byte_range.split("-", 1)[0])
                        if unit != "bytes" or range_start != total:
                            raise RuntimeError(f"unexpected Content-Range: {content_range}")
                        expected_bytes = int(announced_total)
                    elif not total and content_length:
                        expected_bytes = int(content_length)

                    if expected_bytes and expected_bytes > MAX_FILE_BYTES:
                        raise RuntimeError(
                            f"server reports {expected_bytes} bytes, "
                            "above the 100 MB approval limit"
                        )

                    before = total
                    while chunk := response.read(CHUNK_BYTES):
                        total += len(chunk)
                        if total > MAX_FILE_BYTES:
                            raise RuntimeError("download crossed the 100 MB approval limit")
                        digest.update(chunk)
                        stream.write(chunk)

                if expected_bytes is None or total == expected_bytes:
                    break
                if total > expected_bytes:
                    raise RuntimeError(
                        f"response exceeded announced size: received {total} of {expected_bytes}"
                    )
                if total == before:
                    raise RuntimeError("resume attempt made no progress")
            else:
                raise RuntimeError(
                    "response remained truncated after 8 requests: "
                    f"received {total} of {expected_bytes} bytes"
                )

            if expected_bytes is not None and total != expected_bytes:
                raise RuntimeError(
                    f"truncated response: received {total} of {expected_bytes} bytes"
                )

        if not is_pdf(temporary):
            raise RuntimeError("response is not a complete PDF (%PDF- / %%EOF validation failed)")
        temporary.replace(destination)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise

    return digest.hexdigest(), total


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("literature/manifest.csv"),
        help="curated CSV manifest",
    )
    parser.add_argument(
        "--papers-dir",
        type=Path,
        default=Path("literature/papers"),
        help="local ignored PDF directory",
    )
    parser.add_argument("--ids", nargs="*", help="download only these manifest IDs")
    parser.add_argument("--dry-run", action="store_true", help="list selected downloads only")
    parser.add_argument("--overwrite", action="store_true", help="replace existing PDFs")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = load_manifest(args.manifest)
    requested = set(args.ids or [])
    selected = [
        row for row in rows if row["download"] == "1" and (not requested or row["id"] in requested)
    ]
    unknown = requested - {row["id"] for row in rows}
    if unknown:
        print(f"Unknown literature IDs: {', '.join(sorted(unknown))}", file=sys.stderr)
        return 2

    args.papers_dir.mkdir(parents=True, exist_ok=True)
    errors = 0
    for row in selected:
        destination = args.papers_dir / row["local_filename"]
        if args.dry_run:
            print(f"DRY-RUN {row['id']}: {row['pdf_url']} -> {destination}")
            continue

        try:
            if destination.exists() and not args.overwrite:
                if not is_pdf(destination):
                    raise RuntimeError(
                        "existing file is not a complete PDF; rerun with --overwrite"
                    )
                digest, size = sha256_file(destination)
                status = "existing"
            else:
                digest, size = download_pdf(row["pdf_url"], destination)
                status = "downloaded"
            print(
                json.dumps(
                    {
                        "id": row["id"],
                        "status": status,
                        "path": str(destination),
                        "bytes": size,
                        "sha256": digest,
                    },
                    ensure_ascii=False,
                    sort_keys=True,
                )
            )
        except Exception as error:  # noqa: BLE001
            errors += 1
            print(f"ERROR {row['id']}: {error}", file=sys.stderr)

    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
