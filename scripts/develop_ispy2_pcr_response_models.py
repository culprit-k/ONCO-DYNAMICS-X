"""Patient-separated I-SPY 2 response modelling with an RD identifiability gate.

This workflow uses named treatment-stage measurements (T0--T3). It never turns
visit labels into elapsed days or reports day-based kinetic/PDE coefficients.
Only aggregate output is written; participant-level tables stay in memory.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import (
    brier_score_loss,
    log_loss,
    mean_absolute_error,
    mean_squared_error,
    r2_score,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


VOLUME_COLUMNS = {
    "T0": "VOLUME_TUM_BLU_V10",
    "T1": "VOLUME_TUM_BLU_V20",
    "T2": "VOLUME_TUM_BLU_V30",
    "T3": "VOLUME_TUM_BLU_V40",
}
ID_COLUMN = "Patient_ID"
BASE_NUMERIC = ["r1", "r2", "Age_at_Screening", "HR", "HER2"]
QUADRATIC_NUMERIC = BASE_NUMERIC + ["r1_sq", "r1_r2", "r2_sq"]
CATEGORICAL = ["Arm"]


def build_response_frame(features: pd.DataFrame, clinical: pd.DataFrame) -> pd.DataFrame:
    """Strictly link cohorts and derive unit-free log response ratios."""
    feature_key = "CLINICAL-TRIAL-SUBJECT-ID"
    clinical_key = "Patient_ID"
    if feature_key not in features or clinical_key not in clinical:
        raise ValueError("Missing the documented I-SPY 2 participant join key")
    f = features.copy()
    c = clinical.copy()
    f[feature_key] = pd.to_numeric(f[feature_key], errors="raise").astype("Int64")
    c[clinical_key] = pd.to_numeric(c[clinical_key], errors="raise").astype("Int64")
    if f[feature_key].isna().any() or c[clinical_key].isna().any():
        raise ValueError("Participant identifiers cannot be missing")
    if f[feature_key].duplicated().any() or c[clinical_key].duplicated().any():
        raise ValueError("Duplicate participant identifiers; refusing ambiguous linkage")
    linked = f.merge(
        c,
        left_on=feature_key,
        right_on=clinical_key,
        how="left",
        validate="one_to_one",
        indicator=True,
    )
    if not (linked["_merge"] == "both").all():
        raise ValueError("At least one MRI participant lacks a clinical match")
    linked = linked.drop(columns=["_merge"])
    required = [*VOLUME_COLUMNS.values(), "pCR", "Age_at_Screening", "HR", "HER2", "Arm"]
    missing = sorted(set(required) - set(linked.columns))
    if missing:
        raise ValueError(f"Linked data are missing required columns: {missing}")

    linked[ID_COLUMN] = pd.to_numeric(linked[ID_COLUMN], errors="raise").astype("Int64")
    if linked[ID_COLUMN].isna().any() or linked[ID_COLUMN].duplicated().any():
        raise ValueError("The modelling table must contain one complete row per unique patient")

    for stage, column in VOLUME_COLUMNS.items():
        values = pd.to_numeric(linked[column], errors="raise").astype(float)
        if not np.isfinite(values.to_numpy()).all() or (values <= 0).any():
            raise ValueError(f"{column} must contain finite, strictly positive values")
        linked[f"log_volume_{stage}"] = np.log(values)

    linked["r1"] = linked["log_volume_T1"] - linked["log_volume_T0"]
    linked["r2"] = linked["log_volume_T2"] - linked["log_volume_T0"]
    linked["r3"] = linked["log_volume_T3"] - linked["log_volume_T0"]
    linked["r1_sq"] = linked["r1"] ** 2
    linked["r1_r2"] = linked["r1"] * linked["r2"]
    linked["r2_sq"] = linked["r2"] ** 2

    linked["pCR"] = pd.to_numeric(linked["pCR"], errors="raise").astype(int)
    if set(linked["pCR"].unique()) != {0, 1}:
        raise ValueError("pCR must contain both binary values 0 and 1")

    for column in ["Age_at_Screening", "HR", "HER2"]:
        linked[column] = pd.to_numeric(linked[column], errors="raise").astype(float)
        if not np.isfinite(linked[column].to_numpy()).all():
            raise ValueError(f"{column} contains missing or non-finite values")
    if linked["Arm"].isna().any():
        raise ValueError("Treatment arm is missing for at least one linked patient")
    linked["Arm"] = linked["Arm"].astype(str)
    return linked


def _preprocessor(numeric_columns: list[str]) -> ColumnTransformer:
    """Fit all scaling/encoding inside each training-only pipeline."""
    return ColumnTransformer(
        transformers=[
            ("numeric", StandardScaler(), numeric_columns),
            ("arm", OneHotEncoder(handle_unknown="ignore"), CATEGORICAL),
        ],
        remainder="drop",
    )


def _fit_pcr_model(
    train: pd.DataFrame,
    test: pd.DataFrame,
    numeric_columns: list[str],
    c_value: float = 0.5,
) -> tuple[dict[str, float], np.ndarray]:
    columns = numeric_columns + CATEGORICAL
    model = Pipeline(
        steps=[
            ("prepare", _preprocessor(numeric_columns)),
            ("classifier", LogisticRegression(C=c_value, max_iter=5000, solver="lbfgs")),
        ]
    )
    model.fit(train[columns], train["pCR"])
    probability = model.predict_proba(test[columns])[:, 1]
    metrics = {
        "roc_auc": float(roc_auc_score(test["pCR"], probability)),
        "brier_score": float(brier_score_loss(test["pCR"], probability)),
        "log_loss": float(log_loss(test["pCR"], probability)),
    }
    return metrics, probability


def _fit_t3_forecast(
    train: pd.DataFrame,
    test: pd.DataFrame,
    numeric_columns: list[str],
) -> tuple[dict[str, float], np.ndarray]:
    columns = numeric_columns + CATEGORICAL
    model = Pipeline(
        steps=[
            ("prepare", _preprocessor(numeric_columns)),
            ("regressor", Ridge(alpha=5.0)),
        ]
    )
    model.fit(train[columns], train["r3"])
    prediction = model.predict(test[columns])
    metrics = {
        "rmse_log_ratio": float(np.sqrt(mean_squared_error(test["r3"], prediction))),
        "mae_log_ratio": float(mean_absolute_error(test["r3"], prediction)),
        "r2": float(r2_score(test["r3"], prediction)),
    }
    return metrics, prediction


def _bootstrap_comparison(
    y: np.ndarray,
    p_linear: np.ndarray,
    p_quadratic: np.ndarray,
    y_volume: np.ndarray,
    pred_linear: np.ndarray,
    pred_quadratic: np.ndarray,
    *,
    n_bootstrap: int,
    seed: int,
) -> dict[str, Any]:
    if n_bootstrap < 100:
        raise ValueError("n_bootstrap must be at least 100")
    rng = np.random.default_rng(seed)
    auc_linear: list[float] = []
    auc_quadratic: list[float] = []
    delta_auc: list[float] = []
    rmse_linear: list[float] = []
    rmse_quadratic: list[float] = []
    delta_rmse: list[float] = []

    for _ in range(n_bootstrap):
        indices = rng.integers(0, len(y), size=len(y))
        if np.unique(y[indices]).size != 2:
            continue
        auc_l = roc_auc_score(y[indices], p_linear[indices])
        auc_q = roc_auc_score(y[indices], p_quadratic[indices])
        rmse_l = np.sqrt(mean_squared_error(y_volume[indices], pred_linear[indices]))
        rmse_q = np.sqrt(mean_squared_error(y_volume[indices], pred_quadratic[indices]))
        auc_linear.append(float(auc_l))
        auc_quadratic.append(float(auc_q))
        delta_auc.append(float(auc_q - auc_l))
        rmse_linear.append(float(rmse_l))
        rmse_quadratic.append(float(rmse_q))
        delta_rmse.append(float(rmse_q - rmse_l))

    def interval(values: list[float]) -> list[float] | None:
        if not values:
            return None
        return [float(x) for x in np.quantile(values, [0.025, 0.975])]

    return {
        "method": "paired patient bootstrap of the locked test set",
        "requested_resamples": int(n_bootstrap),
        "valid_auc_resamples": int(len(delta_auc)),
        "roc_auc_95ci": {
            "linear": interval(auc_linear),
            "quadratic": interval(auc_quadratic),
        },
        "delta_roc_auc_quadratic_minus_linear_95ci": interval(delta_auc),
        "rmse_log_ratio_95ci": {
            "linear": interval(rmse_linear),
            "quadratic": interval(rmse_quadratic),
        },
        "delta_rmse_quadratic_minus_linear_95ci": interval(delta_rmse),
    }


def reaction_diffusion_identifiability_gate(available: dict[str, bool]) -> dict[str, Any]:
    """Return whether a coupled spatial PDE can be calibrated from the supplied inputs."""
    required = [
        "registered_spatial_tumour_fields",
        "actual_acquisition_times",
        "treatment_exposure_timeline",
        "validated_mri_observation_operator",
        "spatial_microenvironment_fields",
    ]
    missing = [key for key in required if not bool(available.get(key, False))]
    return {
        "status": "eligible_for_spatial_calibration" if not missing else "blocked_nonidentifiable",
        "required_inputs": required,
        "missing_inputs": missing,
        "diffusion_identifiable_from_whole_volume_only": False,
        "mathematical_reason": (
            "For no-flux boundaries, integrating dC/ds = D*laplacian(C) + R over the domain gives "
            "d/ds integral(C) = integral(R), because the boundary integral of grad(C) dot n is zero. "
            "The whole-lesion volume therefore removes D from the aggregate balance; scalar pCR and "
            "morphology summaries do not restore spatial information."
        ),
    }


def analyze(
    features: pd.DataFrame,
    clinical: pd.DataFrame,
    *,
    test_size: float = 0.20,
    seed: int = 20261112,
    n_bootstrap: int = 2000,
) -> dict[str, Any]:
    """Run one locked, patient-level holdout comparison and export aggregates only."""
    if not 0.15 <= test_size <= 0.35:
        raise ValueError("test_size must be between 0.15 and 0.35")
    data = build_response_frame(features, clinical)
    if len(data) < 50:
        raise ValueError("At least 50 linked patients are required for this holdout workflow")

    # Split identifiers (not measurement rows); the wide table has exactly one row/patient.
    patient_ids = data[ID_COLUMN].astype(str).to_numpy()
    labels = data["pCR"].to_numpy()
    train_ids, test_ids = train_test_split(
        patient_ids,
        test_size=test_size,
        random_state=seed,
        stratify=labels,
    )
    train = data[data[ID_COLUMN].astype(str).isin(set(train_ids))].copy()
    test = data[data[ID_COLUMN].astype(str).isin(set(test_ids))].copy()
    if set(train[ID_COLUMN].astype(str)) & set(test[ID_COLUMN].astype(str)):
        raise AssertionError("Patient leakage detected across train/test split")
    if len(train) + len(test) != len(data):
        raise AssertionError("Train/test split did not preserve every linked patient exactly once")

    linear_pcr, p_linear = _fit_pcr_model(train, test, BASE_NUMERIC)
    quadratic_pcr, p_quadratic = _fit_pcr_model(train, test, QUADRATIC_NUMERIC)
    linear_t3, pred_linear = _fit_t3_forecast(train, test, BASE_NUMERIC)
    quadratic_t3, pred_quadratic = _fit_t3_forecast(train, test, QUADRATIC_NUMERIC)

    y_test = test["pCR"].to_numpy(dtype=int)
    y_volume = test["r3"].to_numpy(dtype=float)
    bootstrap = _bootstrap_comparison(
        y_test,
        p_linear,
        p_quadratic,
        y_volume,
        pred_linear,
        pred_quadratic,
        n_bootstrap=n_bootstrap,
        seed=seed,
    )

    available = {
        "registered_spatial_tumour_fields": False,
        "actual_acquisition_times": False,
        "treatment_exposure_timeline": False,
        "validated_mri_observation_operator": False,
        "spatial_microenvironment_fields": False,
    }
    return {
        "analysis": "I-SPY 2 pCR-informed treatment-stage response model comparison",
        "cohort": {
            "feature_rows": int(len(features)),
            "clinical_rows": int(len(clinical)),
            "one_to_one_linked_patients": int(len(data)),
            "feature_unmatched_n": 0,
            "clinical_only_n": int(clinical["Patient_ID"].nunique() - len(data)),
        },
        "outcome": {
            "field": "pCR",
            "positive_code": 1,
            "positive_label": "pCR",
            "negative_code": 0,
            "negative_label": "non-pCR",
            "coding_basis": "Li et al., npj Breast Cancer (2020), DOI 10.1038/s41523-020-00203-7",
            "positive_n": int((data["pCR"] == 1).sum()),
            "negative_n": int((data["pCR"] == 0).sum()),
        },
        "split": {
            "method": "stratified patient-level holdout; split performed on unique participant IDs",
            "test_fraction": float(test_size),
            "seed": int(seed),
            "train_n": int(len(train)),
            "test_n": int(len(test)),
            "train_pCR_positive_n": int((train["pCR"] == 1).sum()),
            "test_pCR_positive_n": int((test["pCR"] == 1).sum()),
            "test_pCR_negative_n": int((test["pCR"] == 0).sum()),
            "patient_overlap_n": 0,
            "all_linked_patients_assigned_once": True,
        },
        "measurement_definition": {
            "primary_measure": "workbook volume field VOLUME_TUM_BLU_V1x; treated as a functional MRI enhancement-volume proxy",
            "physical_unit": "not asserted; workbook-specific unit remains unverified",
            "response_variables": {
                "r1": "ln(V_T1 / V_T0)",
                "r2": "ln(V_T2 / V_T0)",
                "r3": "ln(V_T3 / V_T0)",
            },
            "time_handling": "T1 and T2 are named visit-stage features; no day intervals or per-day rates are constructed",
            "viable_cell_count_equivalence": False,
        },
        "models": {
            "linear_response_reference": {
                "pcr_model": "logit(P(pCR)) = intercept + linear terms in r1, r2, age, HR, HER2 and treatment-arm indicators",
                "t3_forecast": "r3 predicted by ridge regression on r1, r2, age, HR, HER2 and treatment-arm indicators",
                "pcr_metrics_test": linear_pcr,
                "t3_forecast_metrics_test": linear_t3,
            },
            "quadratic_nonlinear_response": {
                "pcr_model": "linear predictor additionally includes r1^2, r1*r2 and r2^2",
                "t3_forecast": "ridge regression additionally includes r1^2, r1*r2 and r2^2",
                "pcr_metrics_test": quadratic_pcr,
                "t3_forecast_metrics_test": quadratic_t3,
            },
        },
        "paired_uncertainty": bootstrap,
        "reaction_diffusion_comparator": reaction_diffusion_identifiability_gate(available),
        "leakage_controls": [
            "One linked row per patient; the split is performed on unique patient identifiers.",
            "No participant identifier is used as a model feature.",
            "The pCR classifier uses T1/T2 only; T3 is not a predictor.",
            "pCR is a training target, never an input feature.",
            "The T3 forecast uses training patients only and does not use pCR as a predictor.",
            "Scaling and treatment-arm encoding are fitted within each training-only pipeline.",
            "Patient-level data and identifiers are not written to the output.",
        ],
        "limitations": [
            "The locked holdout is internal to the 384-person feature-linked cohort, not external validation.",
            "The test set is modest; bootstrap intervals may be wide and are conditional on this split.",
            "The nonlinear term is exploratory; performance differences must be assessed with the paired test-set interval.",
            "Treatment-arm adjustment is predictive, not a causal treatment-effect estimate.",
            "HR and HER2 are used as workbook-coded fields; their detailed assay definitions are not re-derived here.",
            "The workbook volume field's exact unit and implementation equivalence to trial FTV remain unverified.",
            "No spatial PDE coefficient is estimated because the linked workbooks lack the required spatial/time-resolved observations.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", required=True, help="Private local MRI feature workbook")
    parser.add_argument("--clinical", required=True, help="Private local clinical workbook")
    parser.add_argument("--output", required=True, help="Aggregate-only JSON output path")
    parser.add_argument("--test-size", type=float, default=0.20)
    parser.add_argument("--seed", type=int, default=20261112)
    parser.add_argument("--bootstrap", type=int, default=2000)
    args = parser.parse_args()

    features = pd.read_excel(args.features)
    clinical = pd.read_excel(args.clinical)
    result = analyze(
        features,
        clinical,
        test_size=args.test_size,
        seed=args.seed,
        n_bootstrap=args.bootstrap,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "cohort": result["cohort"],
                "split": result["split"],
                "models": result["models"],
                "rd_status": result["reaction_diffusion_comparator"]["status"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
