import csv

from onco_dynamics_x.data_engine.manifest import audit_manifest


def _write_manifest(path, rows):
    columns = [
        "cohort_id", "patient_id", "specimen_id", "sample_id",
        "assay_type", "modality", "data_access", "source_accession",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def test_manifest_audit_passes_valid_rows_and_warns_on_controlled_data(tmp_path):
    path = tmp_path / "manifest.csv"
    _write_manifest(path, [{
        "cohort_id": "cohort-a", "patient_id": "P001", "specimen_id": "T1",
        "sample_id": "S001", "assay_type": "RNA-seq", "modality": "transcriptomics",
        "data_access": "controlled", "source_accession": "ACC-1",
    }])
    report = audit_manifest(path)
    assert report.passed
    assert report.rows == 1
    assert report.unique_patients == 1
    assert report.unique_samples == 1
    assert report.warnings


def test_manifest_audit_detects_duplicate_sample_ids(tmp_path):
    path = tmp_path / "manifest.csv"
    row = {
        "cohort_id": "cohort-a", "patient_id": "P001", "specimen_id": "T1",
        "sample_id": "S001", "assay_type": "RNA-seq", "modality": "transcriptomics",
        "data_access": "open", "source_accession": "ACC-1",
    }
    _write_manifest(path, [row, row])
    report = audit_manifest(path)
    assert not report.passed
    assert any("duplicate sample_id" in error for error in report.errors)


def test_manifest_audit_reports_missing_columns(tmp_path):
    path = tmp_path / "manifest.csv"
    path.write_text("sample_id,patient_id\nS1,P1\n", encoding="utf-8")
    report = audit_manifest(path)
    assert not report.passed
    assert "Missing required columns" in report.errors[0]
