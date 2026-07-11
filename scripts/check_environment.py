"""Report whether the local Python and shared ECG data paths are usable."""

from __future__ import annotations

import importlib.metadata
import platform
import sys

from fo_ekf.paths import resolve_project_paths

REQUIRED_DISTRIBUTIONS = ("numpy", "scipy", "pandas", "matplotlib", "wfdb")


def main() -> int:
    paths = resolve_project_paths()
    print(f"Python: {platform.python_version()} ({sys.executable})")
    print(f"Project root: {paths.root}")
    print(f"ECG data (read-only): {paths.data} [exists={paths.data.exists()}]")
    print(f"Processed data: {paths.processed}")
    print(f"Outputs: {paths.output}")

    missing: list[str] = []
    for distribution in REQUIRED_DISTRIBUTIONS:
        try:
            version = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            missing.append(distribution)
            version = "MISSING"
        print(f"{distribution}: {version}")

    metadata_candidates = (
        paths.data / "raw" / "ptbxl_database.csv",
        paths.data / "ptbxl_database.csv",
    )
    metadata = next((path for path in metadata_candidates if path.exists()), None)
    print(f"PTB-XL metadata: {metadata or 'NOT FOUND'}")

    if missing:
        print(f"ERROR: missing packages: {', '.join(missing)}", file=sys.stderr)
        return 1
    if not paths.data.exists():
        print("ERROR: ECG_DATA_DIR does not exist.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
