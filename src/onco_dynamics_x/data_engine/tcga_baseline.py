"""Paired TCGA-BRCA tumour/normal selection and leakage-safe expression baseline."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from onco_dynamics_x.data_engine.splits import group_train_test_split

PRIMARY_TUMOUR = "Primary Tumor"
SOLID_TISSUE_NORMAL = "Solid Tissue Normal"
_METADATA_FIELDS = [
    "file_id",
    "file_name",
    "md5sum",
    "file_size",
    "access",
    "data_type",
    "workflow_type",
    "case_ids",
    "patient_submitter_ids",
    "sample_submitter_ids",
    "sample_types",
]


def _read_metadata(path: str | Path, expected_type: str) -> dict[str, dict[str, str]]:
    with Path(path).open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        missing = set(_METADATA_FIELDS) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Metadata is missing fields: {', '.join(sorted(missing))}.")
        selected: dict[str, dict[str, str]] = {}
        for row in reader:
            if (row.get("access") or "").lower() != "open":
                raise PermissionError("Only open-access files may enter this workflow.")
            if (row.get("sample_types") or "") != expected_type:
                raise ValueError(
                    f"Expected {expected_type!r} metadata, found {row.get('sample_types')!r}."
                )
            patients = [x for x in (row.get("patient_submitter_ids") or "").split("|") if x]
            samples = [x for x in (row.get("sample_submitter_ids") or "").split("|") if x]
            if len(patients) != 1 or len(samples) != 1:
                raise ValueError(
                    "Each selected expression file must link to exactly one patient and one sample."
                )
            patient_id = patients[0]
            if patient_id not in selected:
                selected[patient_id] = row
    return selected


def select_paired_samples(
    primary_tumour_csv: str | Path,
    normal_csv: str | Path,
    output_metadata_csv: str | Path,
    output_manifest_csv: str | Path,
    max_pairs: int = 20,
) -> list[dict[str, str]]:
    """Select one primary tumour and one normal sample per shared patient.

    The deterministic patient intersection makes each patient a group containing
    both specimen classes. This enables leakage-safe evaluation where a patient's
    tumour and normal samples cannot be split across training and test partitions.
    """
    if max_pairs < 2:
        raise ValueError("max_pairs must be at least 2 for a meaningful grouped evaluation.")
    tumour_by_patient = _read_metadata(primary_tumour_csv, PRIMARY_TUMOUR)
    normal_by_patient = _read_metadata(normal_csv, SOLID_TISSUE_NORMAL)
    paired_patients = sorted(set(tumour_by_patient) & set(normal_by_patient))
    if len(paired_patients) < 8:
        raise ValueError(
            f"Only {len(paired_patients)} paired patients were found; at least 8 are required."
        )
    chosen_patients = paired_patients[:max_pairs]
    metadata_path = Path(output_metadata_csv)
    manifest_path = Path(output_manifest_csv)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    selected_rows: list[dict[str, str]] = []
    manifest_rows: list[dict[str, str]] = []

    for patient_id in chosen_patients:
        for sample_type, source in (
            (PRIMARY_TUMOUR, tumour_by_patient[patient_id]),
            (SOLID_TISSUE_NORMAL, normal_by_patient[patient_id]),
        ):
            row = {field: source.get(field, "") for field in _METADATA_FIELDS}
            selected_rows.append(row)
            manifest_rows.append(
                {
                    "patient_id": patient_id,
                    "sample_id": row["sample_submitter_ids"],
                    "sample_type": sample_type,
                    "label": "1" if sample_type == PRIMARY_TUMOUR else "0",
                    "file_id": row["file_id"],
                    "file_name": row["file_name"],
                    "md5sum": row["md5sum"],
                    "file_size": row["file_size"],
                    "data_access": row["access"],
                }
            )

    with metadata_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_METADATA_FIELDS)
        writer.writeheader()
        writer.writerows(selected_rows)
    manifest_fields = [
        "patient_id", "sample_id", "sample_type", "label",
        "file_id", "file_name", "md5sum", "file_size", "data_access",
    ]
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_fields)
        writer.writeheader()
        writer.writerows(manifest_rows)
    return manifest_rows


def evaluate_tumour_normal_baseline(
    matrix_npz: str | Path,
    selected_manifest_csv: str | Path,
    output_json: str | Path,
    random_seed: int = 42,
    test_size: float = 0.25,
    n_splits: int = 5,
    top_k_genes: int = 100,
) -> dict[str, Any]:
    """Evaluate a modest tumour-vs-normal expression baseline with patient grouping.

    This is a technical pipeline benchmark, not a diagnostic model. All feature
    filtering, selection and scaling are fitted inside each training fold.
    """
    try:
        from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
        from sklearn.linear_model import LogisticRegression
        from sklearn.metrics import (
            balanced_accuracy_score,
            confusion_matrix,
            roc_auc_score,
        )
        from sklearn.model_selection import StratifiedGroupKFold, cross_validate
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler
    except ImportError as error:
        raise RuntimeError(
            "Install the science extra: pip install -e '.[science]'."
        ) from error

    with np.load(matrix_npz, allow_pickle=False) as data:
        sample_ids = tuple(str(x) for x in data["sample_ids"].tolist())
        gene_ids = tuple(str(x) for x in data["gene_ids"].tolist())
        expression = np.asarray(data["values"], dtype=float)
    with Path(selected_manifest_csv).open(
        "r", newline="", encoding="utf-8-sig"
    ) as handle:
        manifest = list(csv.DictReader(handle))
    if not manifest:
        raise ValueError("Selected specimen manifest is empty.")
    by_sample = {row["sample_id"]: row for row in manifest}
    if len(by_sample) != len(manifest):
        raise ValueError("Selected manifest contains duplicate sample identifiers.")
    if set(sample_ids) != set(by_sample):
        raise ValueError("Expression matrix samples do not exactly match the selected manifest.")
    ordered_rows = [by_sample[sample_id] for sample_id in sample_ids]
    patients = np.asarray([row["patient_id"] for row in ordered_rows], dtype=str)
    labels = np.asarray([int(row["label"]) for row in ordered_rows], dtype=int)
    if set(np.unique(labels)) != {0, 1}:
        raise ValueError("Both normal (0) and primary-tumour (1) classes are required.")
    if np.any(expression < 0) or not np.isfinite(expression).all():
        raise ValueError("TPM expression matrix must be finite and non-negative.")
    if len(np.unique(patients)) < n_splits + 1:
        raise ValueError("Not enough unique patients for the requested grouped CV.")
    x = np.log1p(expression)

    train_indices, test_indices = group_train_test_split(
        patients.tolist(), test_size=test_size, random_seed=random_seed
    )
    if set(patients[train_indices]) & set(patients[test_indices]):
        raise RuntimeError("Patient leakage detected between train and test.")
    if set(labels[train_indices]) != {0, 1} or set(labels[test_indices]) != {0, 1}:
        raise ValueError("Grouped holdout must contain both classes in train and test.")

    feature_count = min(top_k_genes, x.shape[1])
    pipeline = Pipeline(
        [
            ("variance", VarianceThreshold(threshold=0.0)),
            ("select", SelectKBest(score_func=f_classif, k=feature_count)),
            ("scale", StandardScaler()),
            (
                "classifier",
                LogisticRegression(
                    class_weight="balanced",
                    max_iter=3000,
                    random_state=random_seed,
                ),
            ),
        ]
    )
    cv = StratifiedGroupKFold(
        n_splits=n_splits, shuffle=True, random_state=random_seed
    )
    cv_results = cross_validate(
        pipeline,
        x,
        labels,
        groups=patients,
        cv=cv,
        scoring={"balanced_accuracy": "balanced_accuracy", "roc_auc": "roc_auc"},
        error_score="raise",
    )
    pipeline.fit(x[train_indices], labels[train_indices])
    test_probabilities = pipeline.predict_proba(x[test_indices])[:, 1]
    test_predictions = (test_probabilities >= 0.5).astype(int)
    matrix = confusion_matrix(labels[test_indices], test_predictions, labels=[0, 1])
    tn, fp, fn, tp = (int(value) for value in matrix.ravel())

    report: dict[str, Any] = {
        "analysis": "paired_primary_tumour_vs_solid_tissue_normal_expression_baseline",
        "cohort": "TCGA-BRCA",
        "assay": "GDC STAR - Counts; tpm_unstranded column",
        "n_samples": len(sample_ids),
        "n_patients": len(np.unique(patients)),
        "class_counts": {
            "solid_tissue_normal": int(np.count_nonzero(labels == 0)),
            "primary_tumour": int(np.count_nonzero(labels == 1)),
        },
        "patient_grouped_holdout": {
            "train_patients": int(len(np.unique(patients[train_indices]))),
            "test_patients": int(len(np.unique(patients[test_indices]))),
            "train_samples": int(len(train_indices)),
            "test_samples": int(len(test_indices)),
            "balanced_accuracy": float(
                balanced_accuracy_score(labels[test_indices], test_predictions)
            ),
            "roc_auc": float(roc_auc_score(labels[test_indices], test_probabilities)),
            "confusion_matrix_labels_normal_then_tumour": {
                "true_normal_predicted_normal": tn,
                "true_normal_predicted_tumour": fp,
                "true_tumour_predicted_normal": fn,
                "true_tumour_predicted_tumour": tp,
            },
            "train_patient_ids": sorted(set(patients[train_indices].tolist())),
            "test_patient_ids": sorted(set(patients[test_indices].tolist())),
        },
        "stratified_group_cv": {
            "folds": n_splits,
            "balanced_accuracy_mean": float(np.mean(cv_results["test_balanced_accuracy"])),
            "balanced_accuracy_std": float(np.std(cv_results["test_balanced_accuracy"], ddof=1)),
            "roc_auc_mean": float(np.mean(cv_results["test_roc_auc"])),
            "roc_auc_std": float(np.std(cv_results["test_roc_auc"], ddof=1)),
            "fold_balanced_accuracy": [
                float(x) for x in cv_results["test_balanced_accuracy"]
            ],
            "fold_roc_auc": [float(x) for x in cv_results["test_roc_auc"]],
        },
        "model": {
            "algorithm": (
                "variance filter + ANOVA top-k + standardization + balanced logistic regression"
            ),
            "top_k_genes": feature_count,
            "input_transform": "log1p(TPM)",
            "feature_selection_fit_within_each_training_fold": True,
            "patient_grouping_applied_to_holdout_and_cross_validation": True,
        },
        "selected_gene_features": [],
        "interpretation": (
            "Technical baseline for distinguishing known specimen labels in this small paired "
            "subset. It is not a screening/diagnostic model, does not predict patient outcomes, "
            "and does not establish treatment response or biological causality."
        ),
    }
    variance_mask = pipeline.named_steps["variance"].get_support()
    retained_genes = np.asarray(gene_ids)[variance_mask]
    selected_mask = pipeline.named_steps["select"].get_support()
    selected_genes = retained_genes[selected_mask]
    coefficients = pipeline.named_steps["classifier"].coef_[0]
    ranked = np.argsort(np.abs(coefficients))[::-1][: min(20, len(selected_genes))]
    report["selected_gene_features"] = [
        {
            "gene_id": str(selected_genes[index]),
            "coefficient": float(coefficients[index]),
        }
        for index in ranked
    ]
    output_path = Path(output_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report
