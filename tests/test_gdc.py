import hashlib
from io import BytesIO
from unittest.mock import patch

import pytest

from onco_dynamics_x.data_engine.gdc import (
    GDCClient,
    GDCFile,
    write_sample_manifest,
)
from onco_dynamics_x.data_engine.manifest import audit_manifest


SAMPLE_FILE_ID = "11111111-2222-3333-4444-555555555555"


def _metadata_response(access="open"):
    return {
        "data": {
            "hits": [
                {
                    "file_id": SAMPLE_FILE_ID,
                    "file_name": "expression.tsv",
                    "md5sum": "d41d8cd98f00b204e9800998ecf8427e",
                    "file_size": 0,
                    "access": access,
                    "data_type": "Gene Expression Quantification",
                    "analysis": {"workflow_type": "STAR - Counts"},
                    "cases": [
                        {
                            "case_id": "case-uuid",
                            "submitter_id": "TCGA-BRCA-01",
                            "samples": [{"submitter_id": "TCGA-BRCA-01-01A"}],
                        }
                    ],
                }
            ],
            "pagination": {"total": 1},
        }
    }


class FakeGDCClient(GDCClient):
    def __init__(self, response):
        super().__init__()
        self.response = response
        self.last_payload = None

    def _post_json(self, endpoint, payload):
        assert endpoint == "files"
        self.last_payload = payload
        return self.response


def test_gdc_open_file_query_uses_public_access_filter_and_parses_linkage():
    client = FakeGDCClient(_metadata_response())
    files = client.query_open_files(workflow_type="STAR - Counts", max_results=10)
    assert len(files) == 1
    assert files[0].access == "open"
    assert files[0].associated_samples == (
        ("case-uuid", "TCGA-BRCA-01", "TCGA-BRCA-01-01A"),
    )
    clauses = client.last_payload["filters"]["content"]
    assert any(
        item["content"]["field"] == "files.access"
        and item["content"]["value"] == "open"
        for item in clauses
    )


def test_gdc_query_rejects_non_open_response():
    client = FakeGDCClient(_metadata_response(access="controlled"))
    with pytest.raises(RuntimeError, match="non-open"):
        client.query_open_files()


def test_gdc_open_file_download_checks_checksum(tmp_path):
    content = b"sample-expression"
    metadata = GDCFile(
        file_id=SAMPLE_FILE_ID,
        file_name="sample.tsv",
        md5sum=hashlib.md5(content, usedforsecurity=False).hexdigest(),
        file_size=len(content),
        access="open",
        data_type="Gene Expression Quantification",
        workflow_type="STAR - Counts",
        case_ids=("case-uuid",),
        patient_submitter_ids=("TCGA-BRCA-01",),
        sample_submitter_ids=("TCGA-BRCA-01-01A",),
        associated_samples=(("case-uuid", "TCGA-BRCA-01", "TCGA-BRCA-01-01A"),),
    )
    with patch("urllib.request.urlopen", return_value=BytesIO(content)):
        path = GDCClient().download_open_file(metadata, tmp_path)
    assert path.read_bytes() == content


def test_gdc_download_rejects_controlled_access(tmp_path):
    metadata = GDCFile(
        file_id=SAMPLE_FILE_ID,
        file_name="sample.tsv",
        md5sum="d41d8cd98f00b204e9800998ecf8427e",
        file_size=0,
        access="controlled",
        data_type="Gene Expression Quantification",
        workflow_type="STAR - Counts",
        case_ids=(),
        patient_submitter_ids=(),
        sample_submitter_ids=(),
        associated_samples=(),
    )
    with pytest.raises(PermissionError):
        GDCClient().download_open_file(metadata, tmp_path)


def test_gdc_manifest_export_passes_sample_audit(tmp_path):
    file = GDCFile(
        file_id=SAMPLE_FILE_ID,
        file_name="sample.tsv",
        md5sum="d41d8cd98f00b204e9800998ecf8427e",
        file_size=0,
        access="open",
        data_type="Gene Expression Quantification",
        workflow_type="STAR - Counts",
        case_ids=("case-uuid",),
        patient_submitter_ids=("TCGA-BRCA-01",),
        sample_submitter_ids=("TCGA-BRCA-01-01A",),
        associated_samples=(("case-uuid", "TCGA-BRCA-01", "TCGA-BRCA-01-01A"),),
    )
    path = write_sample_manifest(
        [file], tmp_path / "gdc_manifest.csv", "TCGA-BRCA", "transcriptomics"
    )
    audit = audit_manifest(path)
    assert audit.passed
    assert audit.unique_patients == 1
    assert audit.unique_samples == 1
