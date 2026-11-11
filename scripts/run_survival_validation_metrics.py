"""Run an age-only TCGA-to-METABRIC survival validation.

The model is fitted in TCGA and frozen before external evaluation. METABRIC
outcomes are used only for evaluation and censoring-weight estimation.
Patient-level records remain in memory and are never written to the report.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from onco_dynamics_x.data_engine.cbioportal import CBioPortalClient
from onco_dynamics_x.data_engine.survival_analysis import fetch_survival_cohort
from onco_dynamics_x.data_engine.survival_metrics import evaluate_ipcw_survival


def _bootstrap_c_index(
    durations: np.ndarray,
    events: np.ndarray,
    risk: np.ndarray,
    *,
    n_bootstrap: int,
    seed: int,
) -> dict[str, Any]:
    from lifelines.utils import concordance_index

    rng = np.random.default_rng(seed)
    estimates: list[float] = []
    n = len(durations)
    for _ in range(n_bootstrap):
        indices = rng.integers(0, n, size=n)
        try:
            value = concordance_index(
                durations[indices], -risk[indices], events[indices]
            )
        except (ZeroDivisionError, ValueError):
            continue
        if np.isfinite(value):
            estimates.append(float(value))
    interval = None
    if len(estimates) >= max(30, int(0.5 * n_bootstrap)):
        interval = [
            float(np.quantile(estimates, 0.025)),
            float(np.quantile(estimates, 0.975)),
        ]
    return {
        "requested_replicates": int(n_bootstrap),
        "valid_replicates": len(estimates),
        "harrell_c_index_95ci": interval,
        "seed": int(seed),
    }


def run_validation(
    *,
    horizons_months: tuple[float, ...] = (12.0, 36.0, 60.0),
    n_bootstrap: int = 250,
    seed: int = 42,
    penalizer: float = 0.1,
) -> dict[str, Any]:
    from lifelines import CoxPHFitter
    from lifelines.statistics import proportional_hazard_test
    from lifelines.utils import concordance_index

    client = CBioPortalClient(timeout=60, retries=3)
    tcga_rows, tcga_qc = fetch_survival_cohort(
        client=client,
        study_id="brca_tcga_pan_can_atlas_2018",
        sample_list_id="brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
        age_attribute="AGE",
        stage_attribute="AJCC_PATHOLOGIC_TUMOR_STAGE",
        stage_level="PATIENT",
        require_stage=False,
    )
    metabric_rows, metabric_qc = fetch_survival_cohort(
        client=client,
        study_id="brca_metabric",
        sample_list_id="brca_metabric_mrna",
        age_attribute="AGE_AT_DIAGNOSIS",
        stage_attribute="TUMOR_STAGE",
        stage_level="SAMPLE",
        require_stage=False,
    )

    tcga = pd.DataFrame(tcga_rows)
    metabric = pd.DataFrame(metabric_rows)
    if tcga.empty or metabric.empty:
        raise RuntimeError("At least one cohort returned no eligible age/OS records.")
    if tcga["patient_id"].duplicated().any() or metabric["patient_id"].duplicated().any():
        raise RuntimeError("Patient grouping failed: duplicate patients remain.")
    if len(tcga) < 100 or int(tcga["event"].sum()) < 30:
        raise ValueError("TCGA development cohort is too small for the prespecified model.")
    if len(metabric) < 100 or int(metabric["event"].sum()) < 30:
        raise ValueError("METABRIC external cohort is too small for evaluation.")
    if not horizons_months or any(value <= 0 for value in horizons_months):
        raise ValueError("Evaluation horizons must be positive.")
    if any(
        right <= left
        for left, right in zip(horizons_months, horizons_months[1:])
    ):
        raise ValueError("Evaluation horizons must be sorted and unique.")

    age_mean = float(tcga["age_years"].mean())
    age_sd = float(tcga["age_years"].std(ddof=0))
    if not np.isfinite(age_sd) or age_sd <= 0:
        raise ValueError("TCGA age standard deviation is invalid.")

    train_covariates = pd.DataFrame(
        {"age_z": (tcga["age_years"] - age_mean) / age_sd}
    )
    test_covariates = pd.DataFrame(
        {"age_z": (metabric["age_years"] - age_mean) / age_sd}
    )
    train = train_covariates.copy()
    train["duration_months"] = tcga["duration_months"].to_numpy()
    train["event"] = tcga["event"].to_numpy()

    model = CoxPHFitter(penalizer=penalizer)
    model.fit(train, duration_col="duration_months", event_col="event")

    test_durations = metabric["duration_months"].to_numpy(dtype=float)
    test_events = metabric["event"].to_numpy(dtype=int)
    test_risk = model.predict_partial_hazard(test_covariates).to_numpy().reshape(-1)
    c_index = float(concordance_index(test_durations, -test_risk, test_events))
    c_index_bootstrap = _bootstrap_c_index(
        test_durations,
        test_events,
        test_risk,
        n_bootstrap=n_bootstrap,
        seed=seed,
    )

    predicted_survival_frame = model.predict_survival_function(
        test_covariates, times=list(horizons_months)
    )
    predicted_survival = predicted_survival_frame.to_numpy(dtype=float).T

    per_horizon_reports = []
    for horizon_index, horizon in enumerate(horizons_months):
        one = evaluate_ipcw_survival(
            test_durations,
            test_events,
            predicted_survival[:, horizon_index],
            horizons_months=(horizon,),
            n_bins=10,
            n_bootstrap=n_bootstrap,
            seed=seed + horizon_index,
        )
        per_horizon_reports.append(one["horizons"][0])

    valid_briers = [
        (float(row["time_months"]), float(row["brier_score_ipcw"]))
        for row in per_horizon_reports
        if row["brier_score_ipcw"] is not None
    ]
    ibs = None
    if len(valid_briers) == len(horizons_months) and len(valid_briers) >= 2:
        grid_times = np.asarray([0.0] + [item[0] for item in valid_briers])
        grid_scores = np.asarray([0.0] + [item[1] for item in valid_briers])
        ibs = float(
            np.sum(
                0.5
                * (grid_scores[1:] + grid_scores[:-1])
                * np.diff(grid_times)
            )
            / grid_times[-1]
        )

    ph_test = proportional_hazard_test(
        model, train, time_transform="rank"
    ).summary
    coefficient = model.summary.loc["age_z"]
    return {
        "analysis": "frozen_age_only_cox_external_survival_validation",
        "endpoint": {
            "name": "overall survival",
            "status_field": "OS_STATUS",
            "time_field": "OS_MONTHS",
            "event_mapping": "1:DECEASED=1; 0:LIVING=0",
            "time_unit": "months",
            "time_origin": (
                "cBioPortal's generic format documents months since initial diagnosis. "
                "Study-specific definitions must be checked against the cohort guides."
            ),
            "nonpositive_times": (
                "Excluded by the cohort assembler and counted in aggregate QC; "
                "no artificial positive time was substituted."
            ),
        },
        "cohorts": {
            "development": {
                "study": "TCGA-BRCA PanCancer Atlas",
                "patients": int(len(tcga)),
                "events": int(tcga["event"].sum()),
                "censored": int((1 - tcga["event"]).sum()),
                "aggregate_qc": tcga_qc,
            },
            "external_evaluation": {
                "study": "METABRIC",
                "patients": int(len(metabric)),
                "events": int(metabric["event"].sum()),
                "censored": int((1 - metabric["event"]).sum()),
                "aggregate_qc": metabric_qc,
            },
        },
        "model": {
            "family": "Cox proportional hazards",
            "predictors": ["age at diagnosis"],
            "stage_used": False,
            "expression_used": False,
            "penalizer": float(penalizer),
            "age_mean_from_tcga": age_mean,
            "age_sd_from_tcga": age_sd,
            "age_hazard_ratio_per_tcga_sd": float(coefficient["exp(coef)"]),
            "age_hazard_ratio_95ci": [
                float(coefficient["exp(coef) lower 95%"]),
                float(coefficient["exp(coef) upper 95%"]),
            ],
            "age_p_value": float(coefficient["p"]),
            "proportional_hazards_rank_test_p_value": float(ph_test.loc["age_z", "p"]),
            "proportional_hazards_test_interpretation": (
                "A non-significant test does not prove proportional hazards."
            ),
        },
        "external_performance": {
            "harrell_c_index": c_index,
            "harrell_c_index_bootstrap": c_index_bootstrap,
            "time_dependent_metrics": {
                "method": (
                    "Cumulative/dynamic IPCW AUC and IPCW Brier; marginal censoring "
                    "Kaplan-Meier estimated in the external evaluation cohort"
                ),
                "censoring_assumption": (
                    "Independent censoring within the external cohort; covariate-dependent "
                    "censoring is not modelled."
                ),
                "horizons": per_horizon_reports,
                "integrated_brier_score_sparse_grid": ibs,
                "integrated_brier_score_note": (
                    "Trapezoidal approximation using time zero and the prespecified "
                    "horizons; not a dense-grid integrated Brier score."
                ),
            },
        },
        "prespecification": {
            "horizons_months": [float(value) for value in horizons_months],
            "calibration_bins": 10,
            "bootstrap_replicates": int(n_bootstrap),
            "bootstrap_seed": int(seed),
            "external_outcomes_used_for_model_fitting_or_tuning": False,
            "stage_mapping_used": False,
        },
        "data_governance": {
            "patient_level_rows_written": False,
            "patient_or_sample_identifiers_written": False,
            "raw_expression_written": False,
            "aggregate_metrics_only": True,
        },
        "interpretation": (
            "Exploratory clinical-only transport evaluation. It is not a clinically "
            "validated prognostic tool, treatment-response model, or validation of the "
            "ONCO-DYNAMICS-X spatial/mechanistic digital twin."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default="results/survival_validation/age_only_external_metrics.json",
    )
    parser.add_argument("--bootstrap", type=int, default=250)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--penalizer", type=float, default=0.1)
    parser.add_argument(
        "--horizons-months",
        type=float,
        nargs="+",
        default=[12.0, 36.0, 60.0],
    )
    args = parser.parse_args()
    if args.bootstrap < 0:
        parser.error("--bootstrap must be non-negative")

    report = run_validation(
        horizons_months=tuple(args.horizons_months),
        n_bootstrap=args.bootstrap,
        seed=args.seed,
        penalizer=args.penalizer,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        "TCGA patients/events:",
        report["cohorts"]["development"]["patients"],
        report["cohorts"]["development"]["events"],
    )
    print(
        "METABRIC patients/events:",
        report["cohorts"]["external_evaluation"]["patients"],
        report["cohorts"]["external_evaluation"]["events"],
    )
    print("METABRIC Harrell C-index:", report["external_performance"]["harrell_c_index"])
    for row in report["external_performance"]["time_dependent_metrics"]["horizons"]:
        print(
            f"{row['time_months']:.0f} months: AUC={row['time_dependent_auc_ipcw']}; "
            f"Brier={row['brier_score_ipcw']}; status={row['metric_status']}"
        )
    print("Aggregate report:", output)


if __name__ == "__main__":
    main()
