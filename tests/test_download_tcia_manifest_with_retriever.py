import json
from pathlib import Path

import pytest

from scripts.download_tcia_manifest_with_retriever import (
    _safe_asset_name,
    run_download,
    validate_uid_manifest,
)


def test_validate_uid_manifest_accepts_unique_series_uids(tmp_path: Path):
    manifest = tmp_path / "manifest.csv"
    manifest.write_text(
        "SeriesInstanceUID\n1.2.840.1\n1.2.840.2\n",
        encoding="utf-8",
    )
    assert validate_uid_manifest(manifest) == ["1.2.840.1", "1.2.840.2"]


def test_validate_uid_manifest_rejects_duplicate_uid(tmp_path: Path):
    manifest = tmp_path / "manifest.csv"
    manifest.write_text(
        "SeriesInstanceUID\n1.2.840.1\n1.2.840.1\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="duplicate"):
        validate_uid_manifest(manifest)


def test_validate_uid_manifest_rejects_wrong_header(tmp_path: Path):
    manifest = tmp_path / "manifest.csv"
    manifest.write_text("series_uid\n1.2.840.1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="header"):
        validate_uid_manifest(manifest)


def test_release_asset_requires_sha256():
    release = {
        "assets": [{
            "name": "TCIA_Data_Retriever_linux_x86_64.zip",
            "browser_download_url": "https://example.invalid/retriever.zip",
        }]
    }
    with pytest.raises(ValueError, match="SHA-256"):
        _safe_asset_name(release)


def test_dry_run_does_not_execute_binary(tmp_path: Path):
    manifest = tmp_path / "manifest.csv"
    manifest.write_text("SeriesInstanceUID\n1.2.840.1\n", encoding="utf-8")
    fake = tmp_path / "retriever"
    fake.write_text("not executed in dry-run", encoding="utf-8")
    fake.chmod(0o755)
    result = run_download(
        manifest,
        tmp_path / "out",
        executable=fake,
        dry_run=True,
    )
    assert result == 0
