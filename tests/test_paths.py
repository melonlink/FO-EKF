from pathlib import Path

from fo_ekf.paths import resolve_project_paths


def test_project_owned_defaults_stay_under_repository(tmp_path: Path) -> None:
    paths = resolve_project_paths(tmp_path)

    assert paths.root == tmp_path.resolve()
    assert paths.processed == (tmp_path / "data" / "processed").resolve()
    assert paths.output == (tmp_path / "artifacts").resolve()


def test_shared_data_default_is_outside_repository(tmp_path: Path) -> None:
    paths = resolve_project_paths(tmp_path)

    assert paths.data == (tmp_path.parent / "PUBLIC" / "data_PTB-XL").resolve()
    assert paths.root not in paths.data.parents
