"""Aggregate-only longitudinal analysis of the I-SPY 2 multi-feature MRI workbook.

The script never exports participant identifiers or row-level measurements.
Visit contrasts are descriptive; no per-day rates or treatment effects are inferred.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats


VOLUME_COLUMNS = [
    "VOLUME_TUM_BLU_V10",
    "VOLUME_TUM_BLU_V20",
    "VOLUME_TUM_BLU_V30",
    "VOLUME_TUM_BLU_V40",
]
CHANGE_COLUMNS = ["FTV_pch_T0_T1", "FTV_pch_T0_T2", "FTV_pch_T0_T3"]
VISITS = ["T0", "T1", "T2", "T3"]
ID_COLUMN = "CLINICAL-TRIAL-SUBJECT-ID"


def _bootstrap_mean(values: np.ndarray, *, iterations: int, rng: np.random.Generator) -> np.ndarray:
    draws = rng.choice(values, size=(iterations, len(values)), replace=True)
    return draws.mean(axis=1)


def analyze_workbook(
    input_path: str | Path, *, bootstrap_iterations: int = 5000, seed: int = 20261112
) -> dict[str, Any]:
    frame = pd.read_excel(input_path)
    required = {ID_COLUMN, *VOLUME_COLUMNS, *CHANGE_COLUMNS}
    missing_columns = required - set(frame.columns)
    if missing_columns:
        raise ValueError(f"Workbook is missing required columns: {sorted(missing_columns)}")

    if frame[ID_COLUMN].isna().any():
        raise ValueError("Missing subject keys found; resolve the source data before analysis.")
    if frame[ID_COLUMN].duplicated().any():
        raise ValueError("Duplicate subject rows found; define a deduplication policy first.")

    volume = frame[VOLUME_COLUMNS].apply(pd.to_numeric, errors="coerce")
    if volume.isna().any().any():
        raise ValueError("Missing/non-numeric volume measurements found.")
    if (volume <= 0).any().any():
        raise ValueError("Non-positive volume values require a censored-data model.")

    # Verify that supplied percent changes are derived from the volume fields.
    formula_errors: dict[str, float] = {}
    for later_index, change_column in zip((1, 2, 3), CHANGE_COLUMNS, strict=True):
        calculated = (volume.iloc[:, later_index] / volume.iloc[:, 0] - 1.0) * 100.0
        error = np.abs(calculated - pd.to_numeric(frame[change_column], errors="coerce"))
        if error.isna().any():
            raise ValueError(f"Missing/non-numeric values in {change_column}.")
        formula_errors[change_column] = float(error.max())
        if formula_errors[change_column] > 1e-6:
            raise ValueError(f"{change_column} does not match the volume-derived percentage change.")

    rng = np.random.default_rng(seed)
    log_volume = np.log(volume.to_numpy(dtype=float))
    n = len(frame)
    visit_summary = []
    for index, visit in enumerate(VISITS):
        values = log_volume[:, index]
        boot = np.exp(_bootstrap_mean(values, iterations=bootstrap_iterations, rng=rng))
        visit_summary.append(
            {
                "visit": visit,
                "n": n,
                "geometric_mean_volume": float(np.exp(values.mean())),
                "geometric_mean_ci95_low": float(np.quantile(boot, 0.025)),
                "geometric_mean_ci95_high": float(np.quantile(boot, 0.975)),
                "median_volume": float(np.median(np.exp(values))),
                "arithmetic_mean_volume": float(np.mean(np.exp(values))),
            }
        )

    contrasts = []
    for index in range(1, 4):
        differences = log_volume[:, index] - log_volume[:, 0]
        boot = _bootstrap_mean(differences, iterations=bootstrap_iterations, rng=rng)
        mean_log_ratio = float(differences.mean())
        contrasts.append(
            {
                "contrast": f"{VISITS[index]}_vs_T0",
                "mean_log_ratio": mean_log_ratio,
                "geometric_mean_ratio": float(np.exp(mean_log_ratio)),
                "ratio_ci95_low": float(np.exp(np.quantile(boot, 0.025))),
                "ratio_ci95_high": float(np.exp(np.quantile(boot, 0.975))),
                "geometric_mean_percent_change": float((np.exp(mean_log_ratio) - 1) * 100),
                "median_individual_percent_change": float(
                    np.median((np.exp(differences) - 1) * 100)
                ),
                "subjects_with_decrease": int((differences < 0).sum()),
                "subjects_with_increase": int((differences > 0).sum()),
            }
        )

    # Reference-prior Normal posterior for subject-level T3/T0 log ratios.
    endpoint = log_volume[:, 3] - log_volume[:, 0]
    posterior_mean = float(endpoint.mean())
    posterior_sd = float(endpoint.std(ddof=1))
    posterior_scale = posterior_sd / np.sqrt(n)
    critical = float(stats.t.ppf(0.975, df=n - 1))
    posterior_low = posterior_mean - critical * posterior_scale
    posterior_high = posterior_mean + critical * posterior_scale

    return {
        "analysis": "aggregate_longitudinal_I_SPY2_MRI_volume",
        "n_subjects": int(n),
        "n_columns": int(frame.shape[1]),
        "missing_cells": int(frame.isna().sum().sum()),
        "duplicate_subject_rows": 0,
        "all_four_volume_values_positive": int((volume > 0).all(axis=1).sum()),
        "source_volume_columns": VOLUME_COLUMNS,
        "percent_change_formula_max_absolute_errors": formula_errors,
        "visit_summary": visit_summary,
        "baseline_contrasts": contrasts,
        "bayesian_endpoint": {
            "likelihood": "Normal subject-level log(T3/T0) contrasts",
            "prior": "reference prior p(mu, sigma^2) proportional to 1/sigma^2",
            "posterior_mean_log_ratio": posterior_mean,
            "posterior_mean_95_credible_interval": [posterior_low, posterior_high],
            "posterior_geometric_ratio": float(np.exp(posterior_mean)),
            "posterior_geometric_ratio_95_credible_interval": [
                float(np.exp(posterior_low)),
                float(np.exp(posterior_high)),
            ],
            "probability_mean_log_ratio_below_zero": float(
                stats.t.cdf(0, df=n - 1, loc=posterior_mean, scale=posterior_scale)
            ),
        },
        "interpretation_limits": [
            "The workbook contains no actual scan acquisition dates or elapsed days.",
            "The workbook contains no treatment-arm/regimen variables.",
            "The workbook contains no HR/HER2 subtype or pCR outcome columns.",
            "The four visits represent unequal treatment intervals; do not infer per-day rates.",
            "Confirm the exact VOLUME_TUM_BLU field definition and physical unit from the source data dictionary.",
            "These are descriptive observational contrasts, not causal treatment effects or clinical predictions.",
        ],
        "privacy": "No participant identifiers or row-level measurements are exported.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Local XLSX workbook.")
    parser.add_argument(
        "--output",
        default="results/ispy2_aggregate_analysis.json",
        help="Aggregate-only JSON output.",
    )
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261112)
    args = parser.parse_args()
    if args.bootstrap < 1000:
        parser.error("--bootstrap must be at least 1000.")

    result = analyze_workbook(args.input, bootstrap_iterations=args.bootstrap, seed=args.seed)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Subjects analysed: {result['n_subjects']}")
    print(f"Missing cells: {result['missing_cells']}")
    print(f"Aggregate report: {output}")


if __name__ == "__main__":
    main()
