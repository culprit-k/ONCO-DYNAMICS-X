"""Censoring-aware external survival metrics for a frozen prediction model.

Implements cumulative/dynamic time-dependent AUC and IPCW Brier scores using
an inverse-probability-of-censoring Kaplan-Meier estimator in the evaluation
cohort. Calibration summaries are descriptive IPCW estimates, not proof of
clinical utility. The marginal censoring estimator assumes independent
censoring within the evaluation cohort.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def _validate_inputs(
    durations: np.ndarray,
    events: np.ndarray,
    predicted_survival: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    durations = np.asarray(durations, dtype=float).reshape(-1)
    events = np.asarray(events, dtype=int).reshape(-1)
    predicted_survival = np.asarray(predicted_survival, dtype=float).reshape(-1)
    if not (len(durations) == len(events) == len(predicted_survival)):
        raise ValueError("durations, events and predictions must have equal lengths.")
    if len(durations) < 2:
        raise ValueError("At least two evaluation patients are required.")
    if not np.all(np.isfinite(durations)) or np.any(durations <= 0):
        raise ValueError("Follow-up times must be finite and strictly positive.")
    if not np.all(np.isin(events, [0, 1])):
        raise ValueError("events must use 1=event and 0=censored.")
    if not np.all(np.isfinite(predicted_survival)):
        raise ValueError("Predicted survival probabilities must be finite.")
    if np.any((predicted_survival < 0) | (predicted_survival > 1)):
        raise ValueError("Predicted survival probabilities must lie in [0, 1].")
    return durations, events, predicted_survival


def _censoring_km(
    durations: np.ndarray,
    events: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return censoring-KM jump times and post-jump survival probabilities."""
    times = np.unique(durations)
    survival = 1.0
    after = []
    for time in times:
        at_risk = int(np.count_nonzero(durations >= time))
        censor_events = int(np.count_nonzero((durations == time) & (events == 0)))
        if at_risk:
            survival *= 1.0 - censor_events / at_risk
        after.append(survival)
    return times, np.asarray(after, dtype=float)


def _g_before(
    km_times: np.ndarray,
    km_after: np.ndarray,
    query_times: np.ndarray,
) -> np.ndarray:
    indices = np.searchsorted(km_times, query_times, side="left") - 1
    values = np.ones(len(query_times), dtype=float)
    valid = indices >= 0
    values[valid] = km_after[indices[valid]]
    return values


def _g_at(
    km_times: np.ndarray,
    km_after: np.ndarray,
    query_time: float,
) -> float:
    index = int(np.searchsorted(km_times, query_time, side="right") - 1)
    return 1.0 if index < 0 else float(km_after[index])


def _weighted_auc(
    case_scores: np.ndarray,
    case_weights: np.ndarray,
    control_scores: np.ndarray,
    control_weights: np.ndarray,
) -> float | None:
    """Weighted Mann-Whitney statistic; ties contribute one half."""
    if len(case_scores) == 0 or len(control_scores) == 0:
        return None
    order = np.argsort(control_scores, kind="mergesort")
    sorted_scores = control_scores[order]
    sorted_weights = control_weights[order]
    cumulative = np.concatenate(([0.0], np.cumsum(sorted_weights)))
    left = np.searchsorted(sorted_scores, case_scores, side="left")
    right = np.searchsorted(sorted_scores, case_scores, side="right")
    concordant_control_weight = cumulative[left] + 0.5 * (
        cumulative[right] - cumulative[left]
    )
    numerator = float(np.sum(case_weights * concordant_control_weight))
    denominator = float(np.sum(case_weights) * np.sum(control_weights))
    if denominator <= 0 or not np.isfinite(denominator):
        return None
    return numerator / denominator


def _one_horizon(
    durations: np.ndarray,
    events: np.ndarray,
    predicted_survival: np.ndarray,
    horizon: float,
    n_bins: int,
    minimum_censoring_survival: float,
) -> dict[str, Any]:
    km_times, km_after = _censoring_km(durations, events)
    g_horizon = _g_at(km_times, km_after, horizon)
    case_mask = (events == 1) & (durations <= horizon)
    control_mask = durations > horizon
    known_mask = case_mask | control_mask
    n_cases = int(case_mask.sum())
    n_controls = int(control_mask.sum())
    n_censored_before = int(len(durations) - known_mask.sum())

    result: dict[str, Any] = {
        "time_months": float(horizon),
        "n_total": int(len(durations)),
        "n_cases": n_cases,
        "n_controls": n_controls,
        "n_censored_before_horizon": n_censored_before,
        "censoring_survival_at_horizon": float(g_horizon),
        "time_dependent_auc_ipcw": None,
        "brier_score_ipcw": None,
        "mean_predicted_survival": float(np.mean(predicted_survival)),
        "ipcw_observed_survival": None,
        "calibration_gap_observed_minus_predicted": None,
        "calibration_bins": [],
        "metric_status": "ok",
    }

    if g_horizon < minimum_censoring_survival:
        result["metric_status"] = (
            "not_estimable: censoring survival at horizon below stability threshold"
        )
        return result
    if n_cases == 0 or n_controls == 0:
        result["metric_status"] = "not_estimable: no cases or no controls at horizon"
        return result

    case_g = _g_before(km_times, km_after, durations[case_mask])
    if np.any(case_g < minimum_censoring_survival):
        result["metric_status"] = (
            "not_estimable: censoring survival before an event below stability threshold"
        )
        return result

    case_weights = 1.0 / case_g
    control_weights = np.full(n_controls, 1.0 / g_horizon, dtype=float)
    risk = 1.0 - predicted_survival
    auc = _weighted_auc(
        risk[case_mask], case_weights, risk[control_mask], control_weights
    )

    known_weights = np.empty(int(known_mask.sum()), dtype=float)
    known_outcomes = np.empty(int(known_mask.sum()), dtype=float)
    known_predictions = predicted_survival[known_mask]
    known_indices = np.flatnonzero(known_mask)
    case_indices = np.flatnonzero(case_mask)
    control_indices = np.flatnonzero(control_mask)
    position = {index: offset for offset, index in enumerate(known_indices)}
    for index, weight in zip(case_indices, case_weights, strict=True):
        offset = position[int(index)]
        known_weights[offset] = weight
        known_outcomes[offset] = 0.0
    for index, weight in zip(control_indices, control_weights, strict=True):
        offset = position[int(index)]
        known_weights[offset] = weight
        known_outcomes[offset] = 1.0

    brier_contributions = known_weights * (known_outcomes - known_predictions) ** 2
    brier = float(np.sum(brier_contributions) / len(durations))
    weight_sum = float(np.sum(known_weights))
    observed_survival = float(np.sum(known_weights * known_outcomes) / weight_sum)
    weighted_prediction = float(np.sum(known_weights * known_predictions) / weight_sum)

    result.update(
        {
            "time_dependent_auc_ipcw": None if auc is None else float(auc),
            "brier_score_ipcw": brier,
            "mean_predicted_survival": weighted_prediction,
            "ipcw_observed_survival": observed_survival,
            "calibration_gap_observed_minus_predicted": (
                observed_survival - weighted_prediction
            ),
        }
    )

    # Equal-sized bins are defined by predicted survival alone, not observed outcomes.
    order = known_indices[np.argsort(predicted_survival[known_indices], kind="mergesort")]
    for bin_indices in np.array_split(order, min(n_bins, len(order))):
        if len(bin_indices) == 0:
            continue
        bin_case_mask = case_mask[bin_indices]
        bin_control_mask = control_mask[bin_indices]
        weights = np.empty(len(bin_indices), dtype=float)
        outcomes = np.empty(len(bin_indices), dtype=float)
        if np.any(bin_case_mask):
            case_rows = bin_indices[bin_case_mask]
            weights[bin_case_mask] = 1.0 / _g_before(
                km_times, km_after, durations[case_rows]
            )
            outcomes[bin_case_mask] = 0.0
        if np.any(bin_control_mask):
            weights[bin_control_mask] = 1.0 / g_horizon
            outcomes[bin_control_mask] = 1.0
        bin_weight_sum = float(np.sum(weights))
        result["calibration_bins"].append(
            {
                "n_evaluable": int(len(bin_indices)),
                "mean_predicted_survival": float(
                    np.sum(weights * predicted_survival[bin_indices]) / bin_weight_sum
                ),
                "ipcw_observed_survival": float(
                    np.sum(weights * outcomes) / bin_weight_sum
                ),
            }
        )
    return result


def evaluate_ipcw_survival(
    durations: np.ndarray,
    events: np.ndarray,
    predicted_survival: np.ndarray,
    horizons_months: list[float] | tuple[float, ...] = (12.0, 36.0, 60.0),
    *,
    n_bins: int = 10,
    n_bootstrap: int = 250,
    seed: int = 42,
    minimum_censoring_survival: float = 0.05,
) -> dict[str, Any]:
    """Evaluate frozen survival predictions at predeclared horizons.

    The censoring distribution is estimated from the evaluation cohort's
    follow-up/censoring indicators for IPCW metric correction. Model fitting
    and prediction must already be complete and must not use evaluation
    outcomes. The censoring KM is marginal and assumes independent censoring.
    """
    durations, events, predicted_survival = _validate_inputs(
        durations, events, predicted_survival
    )
    horizons = np.asarray(horizons_months, dtype=float)
    if horizons.ndim != 1 or len(horizons) == 0:
        raise ValueError("At least one evaluation horizon is required.")
    if not np.all(np.isfinite(horizons)) or np.any(horizons <= 0):
        raise ValueError("Horizons must be finite and strictly positive.")
    if np.any(np.diff(horizons) <= 0):
        raise ValueError("Horizons must be strictly increasing and unique.")
    if n_bins < 2:
        raise ValueError("n_bins must be at least 2.")
    if n_bootstrap < 0:
        raise ValueError("n_bootstrap cannot be negative.")
    if not 0 < minimum_censoring_survival <= 1:
        raise ValueError("minimum_censoring_survival must lie in (0, 1].")

    horizon_results = []
    rng = np.random.default_rng(seed)
    for horizon in horizons:
        point = _one_horizon(
            durations,
            events,
            predicted_survival,
            float(horizon),
            n_bins,
            minimum_censoring_survival,
        )
        auc_samples: list[float] = []
        brier_samples: list[float] = []
        gap_samples: list[float] = []
        if n_bootstrap:
            for _ in range(n_bootstrap):
                indices = rng.integers(0, len(durations), size=len(durations))
                try:
                    boot = _one_horizon(
                        durations[indices],
                        events[indices],
                        predicted_survival[indices],
                        float(horizon),
                        n_bins,
                        minimum_censoring_survival,
                    )
                except (FloatingPointError, ValueError):
                    continue
                if boot["metric_status"] != "ok":
                    continue
                if boot["time_dependent_auc_ipcw"] is not None:
                    auc_samples.append(float(boot["time_dependent_auc_ipcw"]))
                if boot["brier_score_ipcw"] is not None:
                    brier_samples.append(float(boot["brier_score_ipcw"]))
                if boot["calibration_gap_observed_minus_predicted"] is not None:
                    gap_samples.append(
                        float(boot["calibration_gap_observed_minus_predicted"])
                    )

        def interval(values: list[float]) -> list[float] | None:
            if len(values) < max(30, int(0.5 * n_bootstrap)):
                return None
            return [
                float(np.quantile(values, 0.025)),
                float(np.quantile(values, 0.975)),
            ]

        point["bootstrap"] = {
            "requested_replicates": int(n_bootstrap),
            "valid_auc_replicates": len(auc_samples),
            "valid_brier_replicates": len(brier_samples),
            "valid_calibration_gap_replicates": len(gap_samples),
            "time_dependent_auc_ipcw_95ci": interval(auc_samples),
            "brier_score_ipcw_95ci": interval(brier_samples),
            "calibration_gap_95ci": interval(gap_samples),
            "seed": int(seed),
        }
        horizon_results.append(point)

    valid_briers = [
        (float(row["time_months"]), float(row["brier_score_ipcw"]))
        for row in horizon_results
        if row["brier_score_ipcw"] is not None
    ]
    ibs = None
    ibs_interval = None
    if len(valid_briers) == len(horizons) and len(valid_briers) >= 2:
        # Sparse-grid trapezoidal approximation from time zero (Brier=0)
        # through the prespecified evaluation horizons.
        times = np.asarray([0.0] + [item[0] for item in valid_briers])
        scores = np.asarray([0.0] + [item[1] for item in valid_briers])
        ibs = float(np.trapezoid(scores, times) / times[-1])

    return {
        "method": (
            "Cumulative/dynamic IPCW AUC and IPCW Brier score; marginal "
            "Kaplan-Meier censoring model estimated in the external evaluation cohort"
        ),
        "assumptions": {
            "event_encoding": "1=death/event; 0=right-censored",
            "time_unit": "months",
            "censoring": "marginal independent censoring within the evaluation cohort",
            "evaluation_outcomes_used_for_model_fitting": False,
            "external_censoring_distribution_used_for_ipcw_only": True,
            "minimum_censoring_survival": float(minimum_censoring_survival),
        },
        "n_evaluation_patients": int(len(durations)),
        "integrated_brier_score_sparse_grid": ibs,
        "integrated_brier_score_interval": ibs_interval,
        "integrated_brier_score_note": (
            "Trapezoidal approximation over time zero and the prespecified horizons; "
            "not a dense-grid integrated Brier score."
        ),
        "horizons": horizon_results,
    }
