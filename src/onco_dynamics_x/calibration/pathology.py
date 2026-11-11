"""Cross-validated association of expression modules with pathology annotations.

This evaluates whether descriptive molecular modules track available
cross-sectional pathology labels. It does not estimate PDE coefficients,
validate an imaging observation operator, or establish causal mechanisms.
"""
from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


def evaluate_module_pathology_association(
    features: NDArray[np.float64],
    patient_ids: list[str] | tuple[str, ...],
    labels_by_patient: dict[str, str],
    *,
    module_names: tuple[str, ...] | list[str],
    n_splits: int = 5,
    n_bootstrap: int = 500,
    seed: int = 42,
) -> dict[str, Any]:
    """Generate patient-level out-of-fold classification metrics.

    The fixed model is StandardScaler + class-balanced L2 logistic regression
    (C=1). There is no hyperparameter tuning. Patient IDs are used only for
    alignment/deduplication and are never included in the returned report.
    """
    try:
        from sklearn.linear_model import LogisticRegression
        from sklearn.metrics import (
            accuracy_score,
            balanced_accuracy_score,
            confusion_matrix,
            f1_score,
            log_loss,
            precision_recall_fscore_support,
        )
        from sklearn.model_selection import StratifiedKFold, cross_val_predict
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler
    except ImportError as error:
        raise ImportError("Install the science extra: pip install -e '.[science]'.") from error

    matrix = np.asarray(features, dtype=float)
    ids = tuple(str(value) for value in patient_ids)
    names = tuple(str(value) for value in module_names)
    if matrix.ndim != 2 or matrix.shape != (len(ids), len(names)):
        raise ValueError("features must align with patient_ids and module_names.")
    if not np.isfinite(matrix).all():
        raise ValueError("features must contain only finite values.")
    if len(set(ids)) != len(ids):
        raise ValueError("patient_ids must be unique; aggregate repeated samples first.")
    if len(set(names)) != len(names) or not names:
        raise ValueError("module_names must be non-empty and unique.")
    if n_splits < 2 or n_bootstrap < 0:
        raise ValueError("n_splits must be at least two and n_bootstrap non-negative.")

    selected = [
        index for index, patient_id in enumerate(ids)
        if patient_id in labels_by_patient
        and str(labels_by_patient[patient_id]).strip()
        and str(labels_by_patient[patient_id]).strip().upper()
        not in {"NA", "N/A", "UNKNOWN", "UNDEFINED", "UNDEF", "NOT REPORTED"}
    ]
    if len(selected) < n_splits * 2:
        raise ValueError("Insufficient labelled patients for cross-validation.")
    x = matrix[selected]
    y = np.asarray([str(labels_by_patient[ids[index]]).strip() for index in selected])
    classes, class_counts = np.unique(y, return_counts=True)
    if len(classes) < 2 or int(class_counts.min()) < n_splits:
        raise ValueError("Each pathology class must contain at least n_splits patients.")

    estimator = make_pipeline(
        StandardScaler(),
        LogisticRegression(
            C=1.0,
            class_weight="balanced",
            solver="lbfgs",
            max_iter=3000,
        ),
    )
    splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    predicted = cross_val_predict(estimator, x, y, cv=splitter, method="predict")
    probabilities = cross_val_predict(
        estimator, x, y, cv=splitter, method="predict_proba"
    )
    class_order = np.unique(y).tolist()
    one_hot = np.column_stack([y == label for label in class_order]).astype(float)
    brier = float(np.mean(np.sum((probabilities - one_hot) ** 2, axis=1)))
    accuracy = float(accuracy_score(y, predicted))
    balanced = float(balanced_accuracy_score(y, predicted))
    macro_f1 = float(f1_score(y, predicted, average="macro", zero_division=0))
    logloss = float(log_loss(y, probabilities, labels=class_order))
    precision, recall, f1, support = precision_recall_fscore_support(
        y, predicted, labels=class_order, zero_division=0
    )
    majority_accuracy = float(class_counts.max() / class_counts.sum())

    rng = np.random.default_rng(seed)
    class_indices = [np.flatnonzero(y == label) for label in class_order]
    bootstrap_metrics = {"accuracy": [], "balanced_accuracy": [], "macro_f1": []}
    for _ in range(n_bootstrap):
        sample_indices = np.concatenate(
            [
                rng.choice(indices, size=len(indices), replace=True)
                for indices in class_indices
            ]
        )
        y_boot = y[sample_indices]
        pred_boot = predicted[sample_indices]
        bootstrap_metrics["accuracy"].append(
            float(accuracy_score(y_boot, pred_boot))
        )
        bootstrap_metrics["balanced_accuracy"].append(
            float(balanced_accuracy_score(y_boot, pred_boot))
        )
        bootstrap_metrics["macro_f1"].append(
            float(f1_score(y_boot, pred_boot, average="macro", zero_division=0))
        )
    intervals = {
        metric: (
            [
                float(np.quantile(values, 0.025)),
                float(np.quantile(values, 0.975)),
            ]
            if values
            else None
        )
        for metric, values in bootstrap_metrics.items()
    }

    return {
        "patients_evaluated": int(len(selected)),
        "patients_excluded_missing_or_unknown_label": int(len(ids) - len(selected)),
        "module_names": list(names),
        "classes": class_order,
        "class_counts": {
            label: int(count) for label, count in zip(class_order, class_counts, strict=True)
        },
        "cross_validation": {
            "method": f"{n_splits}-fold stratified patient-level out-of-fold predictions",
            "hyperparameter_tuning": False,
            "model": "StandardScaler + class-balanced L2 logistic regression (C=1)",
        },
        "metrics": {
            "accuracy": accuracy,
            "balanced_accuracy": balanced,
            "macro_f1": macro_f1,
            "multiclass_brier_score": brier,
            "log_loss": logloss,
            "bootstrap_95ci": intervals,
        },
        "majority_class_baseline_accuracy": majority_accuracy,
        "per_class": {
            label: {
                "precision": float(precision[index]),
                "recall": float(recall[index]),
                "f1": float(f1[index]),
                "support": int(support[index]),
            }
            for index, label in enumerate(class_order)
        },
        "confusion_matrix": {
            "rows_actual_columns_predicted": confusion_matrix(
                y, predicted, labels=class_order
            ).tolist(),
            "class_order": class_order,
        },
        "interpretation": (
            "cross-sectional association with existing cohort pathology labels; "
            "not a causal mechanism, PDE parameter estimate, or clinical predictor"
        ),
        "data_governance": {
            "patient_ids_exported": False,
            "row_level_labels_exported": False,
            "out_of_fold_predictions_exported": False,
        },
    }
