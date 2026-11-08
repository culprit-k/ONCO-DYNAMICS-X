import csv

import numpy as np
import pytest

from onco_dynamics_x.data_engine.tcga_baseline import (
    evaluate_tumour_normal_baseline,
    select_paired_samples,
)

FIELDS = [
    "file_id", "file_name", "md5sum", "file_size", "access", "data_type",
    "workflow_type", "case_ids", "patient_submitter_ids",
    "sample_submitter_ids", "sample_types",
]


def _write_rows(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def _row(patient, sample_type, suffix):
    return {
        "file_id": f"file-{patient}-{suffix}",
        "file_name": f"{patient}-{suffix}.tsv",
        "md5sum": "a" * 32,
        "file_size": "10",
        "access": "open",
        "data_type": "Gene Expression Quantification",
        "workflow_type": "STAR - Counts",
        "case_ids": f"case-{patient}",
        "patient_submitter_ids": patient,
        "sample_submitter_ids": f"{patient}-{suffix}",
        "sample_types": sample_type,
    }


def test_select_paired_samples_keeps_one_sample_per_patient_and_class(tmp_path):
    tumour = tmp_path / "tumour.csv"
    normal = tmp_path / "normal.csv"
    _write_rows(tumour, [_row(f"P{i}", "Primary Tumor", "T") for i in range(10)])
    _write_rows(normal, [_row(f"P{i}", "Solid Tissue Normal", "N") for i in range(2, 12)])
    metadata = tmp_path / "selected.csv"
    manifest = tmp_path / "manifest.csv"
    rows = select_paired_samples(tumour, normal, metadata, manifest, max_pairs=5)
    assert len(rows) == 10
    assert len({row["patient_id"] for row in rows}) == 5
    assert {row["sample_type"] for row in rows} == {
        "Primary Tumor", "Solid Tissue Normal"
    }
    with manifest.open(newline="", encoding="utf-8") as handle:
        saved = list(csv.DictReader(handle))
    assert len(saved) == 10
    assert all(row["data_access"] == "open" for row in saved)


def test_select_paired_samples_rejects_non_open_access(tmp_path):
    tumour = tmp_path / "tumour.csv"
    normal = tmp_path / "normal.csv"
    rows_t = [_row(f"P{i}", "Primary Tumor", "T") for i in range(10)]
    rows_n = [_row(f"P{i}", "Solid Tissue Normal", "N") for i in range(10)]
    rows_t[0]["access"] = "controlled"
    _write_rows(tumour, rows_t)
    _write_rows(normal, rows_n)
    with pytest.raises(PermissionError):
        select_paired_samples(
            tumour, normal, tmp_path / "selected.csv", tmp_path / "manifest.csv"
        )


def test_tumour_normal_evaluation_groups_patients_and_writes_report(tmp_path):
    pytest.importorskip("sklearn")
    samples = []
    genes = [f"ENSG{i:05d}" for i in range(60)]
    values = []
    for i in range(20):
        patient = f"P{i:03d}"
        normal_id = f"{patient}-N"
        tumour_id = f"{patient}-T"
        samples.extend(
            [
                {
                    "patient_id": patient,
                    "sample_id": normal_id,
                    "sample_type": "Solid Tissue Normal",
                    "label": "0",
                },
                {
                    "patient_id": patient,
                    "sample_id": tumour_id,
                    "sample_type": "Primary Tumor",
                    "label": "1",
                },
            ]
        )
        rng = np.random.default_rng(i)
        normal = rng.lognormal(mean=0.5, sigma=0.2, size=len(genes))
        tumour = rng.lognormal(mean=0.5, sigma=0.2, size=len(genes))
        tumour[:8] *= 8.0
        values.extend([normal, tumour])
    matrix_path = tmp_path / "matrix.npz"
    np.savez_compressed(
        matrix_path,
        sample_ids=np.asarray([row["sample_id"] for row in samples]),
        gene_ids=np.asarray(genes),
        values=np.asarray(values),
    )
    manifest_path = tmp_path / "manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(samples[0]))
        writer.writeheader()
        writer.writerows(samples)
    output = tmp_path / "report.json"
    report = evaluate_tumour_normal_baseline(
        matrix_path, manifest_path, output, n_splits=5, top_k_genes=20
    )
    assert report["n_patients"] == 20
    assert report["n_samples"] == 40
    assert report["patient_grouped_holdout"]["train_patients"] == 15
    assert report["patient_grouped_holdout"]["test_patients"] == 5
    assert set(report["patient_grouped_holdout"]["train_patient_ids"]).isdisjoint(
        report["patient_grouped_holdout"]["test_patient_ids"]
    )
    assert output.exists()
