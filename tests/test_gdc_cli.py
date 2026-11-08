import csv
from argparse import Namespace

import pytest

from onco_dynamics_x import cli
from onco_dynamics_x.data_engine.gdc import GDCFile


def _write_metadata(path, access="open"):
    columns = [
        "file_id", "file_name", "md5sum", "file_size", "access",
        "data_type", "workflow_type", "case_ids",
        "patient_submitter_ids", "sample_submitter_ids",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerow({
            "file_id": "11111111-2222-3333-4444-555555555555",
            "file_name": "expression.tsv",
            "md5sum": "d41d8cd98f00b204e9800998ecf8427e",
            "file_size": 0,
            "access": access,
            "data_type": "Gene Expression Quantification",
            "workflow_type": "STAR - Counts",
            "case_ids": "case-1",
            "patient_submitter_ids": "patient-1",
            "sample_submitter_ids": "sample-1",
        })


def test_gdc_download_command_refuses_controlled_metadata(tmp_path):
    metadata = tmp_path / "files.csv"
    _write_metadata(metadata, access="controlled")
    args = Namespace(
        metadata_csv=str(metadata),
        file_id="11111111-2222-3333-4444-555555555555",
        destination=str(tmp_path / "raw"),
        timeout=5.0,
    )
    with pytest.raises(PermissionError, match="open-access"):
        cli.run_gdc_download_open(args)


def test_gdc_download_command_uses_verified_metadata_row(tmp_path, monkeypatch):
    metadata = tmp_path / "files.csv"
    _write_metadata(metadata)
    downloaded_path = tmp_path / "raw" / "expression.tsv"

    class FakeClient:
        def __init__(self, timeout):
            assert timeout == 5.0

        def download_open_file(self, record, destination):
            assert isinstance(record, GDCFile)
            assert record.access == "open"
            assert record.sample_submitter_ids == ("sample-1",)
            assert destination == str(tmp_path / "raw")
            return downloaded_path

    monkeypatch.setattr(cli, "GDCClient", FakeClient)
    args = Namespace(
        metadata_csv=str(metadata),
        file_id="11111111-2222-3333-4444-555555555555",
        destination=str(tmp_path / "raw"),
        timeout=5.0,
    )
    assert cli.run_gdc_download_open(args) == downloaded_path
