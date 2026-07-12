"""Run the locked four-record Fantasia feasibility pilot."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from fo_ekf.fantasia_pilot import FantasiaPilotConfig, run_fantasia_pilot


def _default_data_dir() -> Path:
    raw = os.environ.get("ECG_DATA_DIR")
    if raw:
        return Path(raw)
    return Path(__file__).resolve().parents[2] / "PUBLIC" / "fantasia"


def _default_output_dir() -> Path:
    raw = os.environ.get("ECG_OUTPUT_DIR")
    if raw:
        return Path(raw) / "fantasia_pilot_candidate"
    return Path(__file__).resolve().parents[1] / "artifacts" / "fantasia_pilot_candidate"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=_default_data_dir())
    parser.add_argument("--output-dir", type=Path, default=_default_output_dir())
    arguments = parser.parse_args()
    summary = run_fantasia_pilot(
        arguments.data_dir,
        arguments.output_dir,
        config=FantasiaPilotConfig(),
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
