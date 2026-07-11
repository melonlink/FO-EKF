"""Resolve project paths without mutating the shared source-data directory."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ProjectPaths:
    """Input and output locations used by experiments."""

    root: Path
    data: Path
    processed: Path
    output: Path


def resolve_project_paths(root: Path | None = None) -> ProjectPaths:
    """Resolve paths from environment variables and repository-safe defaults."""

    project_root = (root or Path(__file__).resolve().parents[2]).resolve()

    def resolve_value(variable: str, default: Path) -> Path:
        raw_value = os.environ.get(variable)
        value = Path(raw_value) if raw_value else default
        if not value.is_absolute():
            value = project_root / value
        return value.resolve()

    return ProjectPaths(
        root=project_root,
        data=resolve_value("ECG_DATA_DIR", project_root.parent / "PUBLIC" / "data_PTB-XL"),
        processed=resolve_value("ECG_PROCESSED_DIR", project_root / "data" / "processed"),
        output=resolve_value("ECG_OUTPUT_DIR", project_root / "artifacts"),
    )
