"""Aggregate-only longitudinal analysis of the I-SPY 2 multi-feature MRI workbook.

This analysis is descriptive/exploratory. It does not treat nominal MRI visits as
elapsed time, infer causal treatment effects, or estimate mechanistic PDE rates.
No participant identifiers or row-level observations are exported.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

ID_COLUMN = "CLINICAL-TRIAL-SUBJECT-ID"
FEATURES: dict[str, tuple[str, ...]] = {
    "volume_metric_candidate_ftv": (
        "VOLUME_TUM_BLU_V10",
        "VOLUME_TUM_BLU_V20",
        "VOLUME_TUM_BLU_V30",
        "VOLUME_TUM_BLU_V40",
    ),
    "sphericity": ("SPHERICITY_T0", "SPHERICITY_T1", "SPHERICITY_T2", "SPHERICITY_T3"),
    "longest_diameter": ("LD_T0", "LD_T1", "LD_T2", "LD_T3"),
    "background_parenchymal_enhancement": (
        "BPE_5slice_mean_T0",
        "BPE_5slice_mean_T1",
        "BPE_5slice_mean_T2",
        "BPE_5slice_mean_T3",
    ),
}
CHANGE_COLUMNS: dict[str, tuple[str, ...]] = {
    "volume_metric_candidate_ftv": (
        "FTV_pch_T0_T1", "FTV_pch_T0_T2", "FTV_pch_T0_T3"
    ),
    "sphericity": (
        "Sphericity_pch_T0_T1", "Sphericity_pch_T0_T2", "Sphericity_pch_T0_T3"
    ),
    "longest_diameter": ("LD_pch_T0_T1", "LD_pch_T0_T2", "LD_pch_T0_T3"),
    "background_parenchymal_enhancement": (
        "BPE_pch_T0_T1", "BPE_pch_T0_T2", "BPE_pch_T0_T3"
    ),
}
VISITS = ("T0", "T1", "T2", "T3")


def _holm_adjust(p_values: list[float]) -> list[float]:
    """Holm family-wise adjusted p-values, returned in original order."""
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    running = 0.0
    m = len(p)
    for rank, idx in enumerate(order):
        candidate = (m - rank) * p[idx]
        running = max(running, candidate)
        adjusted[idx] = min(1.0, running)
    return adjusted.tolist()


def _bootstrap_median_ci(
    values: np.ndarray,
    *,
    iterations: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    if len(values) == 0:
        return float("nan"), float("nan")
    draws = rng.choice(values, size=(iterations, len(values)), replace=True)
    medians = np.median(draws, axis=1)
    low, high = np.quantile(medians, [0.025, 0.975])
    return float(low), float(high)


def _bootstrap_spearman_ci(
    x: np.ndarray,
    y: np.ndarray,
    *,
    iterations: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    estimates: list[float] = []
    n = len(x)
    for _ in range(iterations):
        idx = rng.integers(0, n, size=n)
        if np.ptp(x[idx]) == 0 or np.ptp(y[idx]) == 0:
            continue
        rho = stats.spearmanr(x[idx], y[idx]).statistic
        if np.isfinite(rho):
            estimates.append(float(rho))
    if not estimates:
        return float("nan"), float("nan")
    low, high = np.quantile(estimates, [0.025, 0.975])
    return float(low), float(high)


def analyze_frame(
    frame: pd.DataFrame,
    *,
    bootstrap_iterations: int = 2000,
    seed: int = 20261112,
) -> dict[str, Any]:
    if bootstrap_iterations < 200:
        raise ValueError("Use at least 200 bootstrap iterations.")
    required = {ID_COLUMN}
    for cols in FEATURES.values():
        required.update(cols)
    for cols in CHANGE_COLUMNS.values():
        required.update(cols)
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Workbook is missing required columns: {sorted(missing)}")
    if frame[ID_COLUMN].isna().any() or frame[ID_COLUMN].duplicated().any():
        raise ValueError("Subject identifiers must be present and unique for paired analysis.")
    if frame.empty:
        raise ValueError("No participant records found.")

    rng = np.random.default_rng(seed)
    raw: dict[str, np.ndarray] = {}
    formula_errors: dict[str, float] = {}
    for feature, cols in FEATURES.items():
        matrix = frame[list(cols)].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        if not np.isfinite(matrix).all():
            raise ValueError(f"Non-finite or missing values in feature {feature}.")
        if (matrix < 0).any():
            raise ValueError(f"Negative measurements found in feature {feature}.")
        if (matrix[:, 0] <= 0).any():
            raise ValueError(f"Baseline measurements must be positive for {feature}.")
        raw[feature] = matrix
        for visit_index, change_col in enumerate(CHANGE_COLUMNS[feature], start=1):
            observed = pd.to_numeric(frame[change_col], errors="coerce").to_numpy(float)
            calculated = (matrix[:, visit_index] / matrix[:, 0] - 1.0) * 100.0
            if not np.isfinite(observed).all():
                raise ValueError(f"Non-finite or missing values in {change_col}.")
            max_error = float(np.max(np.abs(calculated - observed)))
            formula_errors[change_col] = max_error
            if max_error > 1e-6:
                raise ValueError(f"{change_col} fails the derived-percentage formula check.")

    n = len(frame)
    feature_summaries: dict[str, Any] = {}
    global_friedman_p: list[float] = []
    posthoc_rows: list[dict[str, Any]] = []

    for feature, matrix in raw.items():
        visit_summary = []
        for j, visit in enumerate(VISITS):
            values = matrix[:, j]
            q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
            item: dict[str, Any] = {
                "visit": visit,
                "n": n,
                "mean": float(np.mean(values)),
                "median": float(median),
                "q1": float(q1),
                "q3": float(q3),
                "minimum": float(np.min(values)),
                "maximum": float(np.max(values)),
                "zero_count": int(np.sum(values == 0)),
            }
            if feature == "volume_metric_candidate_ftv":
                log_values = np.log(values)
                boot = np.exp(
                    rng.choice(
                        log_values, size=(bootstrap_iterations, n), replace=True
                    ).mean(axis=1)
                )
                item["geometric_mean"] = float(np.exp(log_values.mean()))
                item["geometric_mean_ci95"] = [
                    float(v) for v in np.quantile(boot, [0.025, 0.975])
                ]
            visit_summary.append(item)

        friedman = stats.friedmanchisquare(*[matrix[:, j] for j in range(4)])
        global_friedman_p.append(float(friedman.pvalue))
        feature_summaries[feature] = {
            "visits": visit_summary,
            "friedman_statistic": float(friedman.statistic),
            "friedman_p_value": float(friedman.pvalue),
            "interpretation": "overall within-subject difference across four ordered visit stages",
        }

        for j in range(1, 4):
            before = matrix[:, 0]
            after = matrix[:, j]
            pct = (after / before - 1.0) * 100.0
            delta = after - before
            test = stats.wilcoxon(after, before, zero_method="zsplit", alternative="two-sided")
            pct_ci = _bootstrap_median_ci(pct, iterations=bootstrap_iterations, rng=rng)
            delta_ci = _bootstrap_median_ci(delta, iterations=bootstrap_iterations, rng=rng)
            posthoc_rows.append(
                {
                    "feature": feature,
                    "contrast": f"{VISITS[j]}_vs_T0",
                    "median_absolute_change": float(np.median(delta)),
                    "median_absolute_change_ci95": list(delta_ci),
                    "median_percent_change": float(np.median(pct)),
                    "median_percent_change_ci95": list(pct_ci),
                    "n_decreased": int(np.sum(delta < 0)),
                    "n_increased": int(np.sum(delta > 0)),
                    "n_unchanged": int(np.sum(delta == 0)),
                    "wilcoxon_statistic": float(test.statistic),
                    "p_value": float(test.pvalue),
                }
            )

    omnibus_adjusted = _holm_adjust(global_friedman_p)
    for feature, adjusted in zip(FEATURES, omnibus_adjusted, strict=True):
        feature_summaries[feature]["friedman_p_holm_4_features"] = adjusted
    posthoc_adjusted = _holm_adjust([row["p_value"] for row in posthoc_rows])
    for row, adjusted in zip(posthoc_rows, posthoc_adjusted, strict=True):
        row["p_value_holm_12_contrasts"] = adjusted

    volume_pct = pd.to_numeric(frame["FTV_pch_T0_T3"], errors="coerce").to_numpy(float)
    associations: list[dict[str, Any]] = []
    for feature, change_col in (
        ("sphericity", "Sphericity_pch_T0_T3"),
        ("longest_diameter", "LD_pch_T0_T3"),
        ("background_parenchymal_enhancement", "BPE_pch_T0_T3"),
    ):
        other_pct = pd.to_numeric(frame[change_col], errors="coerce").to_numpy(float)
        valid = np.isfinite(volume_pct) & np.isfinite(other_pct)
        rho, p_value = stats.spearmanr(volume_pct[valid], other_pct[valid])
        ci = _bootstrap_spearman_ci(
            volume_pct[valid], other_pct[valid], iterations=bootstrap_iterations, rng=rng
        )
        associations.append(
            {
                "feature_change": feature,
                "paired_n": int(valid.sum()),
                "spearman_rho": float(rho),
                "spearman_ci95_bootstrap": list(ci),
                "p_value_exploratory": float(p_value),
            }
        )

    return {
        "analysis": "ispy2_multifeature_longitudinal_aggregate_analysis",
        "n_participants": n,
        "visits": list(VISITS),
        "privacy": "aggregate-only; no subject identifiers or row-level observations exported",
        "measurement_mapping_caveat": (
            "VOLUME_TUM_BLU_V10-V40 are treated as the workbook's longitudinal volume field. "
            "The workbook headers alone do not establish physical units or prove that this field "
            "is interchangeable with a DICOM-derived FTV implementation. Verify the workbook "
            "specific data dictionary before assigning units or a physical PDE observation operator."
        ),
        "temporal_caveat": (
            "T0-T3 are treatment-stage visits with unequal intervals. No per-day rates, "
            "untreated growth rates, causal treatment effects, or drug-response parameters are estimated."
        ),
        "formula_max_absolute_error_percentage_points": formula_errors,
        "feature_summaries": feature_summaries,
        "baseline_contrasts": posthoc_rows,
        "t3_change_associations_with_volume_field": associations,
        "multiple_testing": {
            "omnibus": "Holm adjustment across four feature-level Friedman tests",
            "baseline_contrasts": "Holm adjustment across twelve feature-by-visit Wilcoxon tests",
            "correlations": (
                "exploratory; three associations, report nominal p-values and bootstrap CIs"
            ),
        },
        "limits": [
            "Visit-stage comparisons are not elapsed-time kinetic estimates.",
            "All participants are treatment-exposed; no untreated tumour growth rate "
            "is identified.",
            "The volume-field physical unit and exact mapping require the "
            "workbook-specific dictionary.",
            "FTV-like enhancement volume is not direct viable tumour-cell burden.",
            "Associations among feature changes are exploratory and not causal or mechanistic.",
            "No pCR, subtype, treatment arm, actual acquisition dates or drug "
            "exposure data are present in this workbook.",
            "No model parameter in the ONCO-DYNAMICS-X PDE is estimated by this analysis.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Local workbook path; remains private.")
    parser.add_argument(
        "--output",
        default="results/ispy2_multifeature_aggregate.json",
        help="Aggregate JSON output.",
    )
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261112)
    args = parser.parse_args()
    if args.bootstrap < 200:
        parser.error("--bootstrap must be at least 200.")
    frame = pd.read_excel(args.input)
    report = analyze_frame(frame, bootstrap_iterations=args.bootstrap, seed=args.seed)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Participants analysed: {report['n_participants']}")
    print(f"Aggregate output: {output_path}")


if __name__ == "__main__":
    main()
