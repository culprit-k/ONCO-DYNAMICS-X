import numpy as np
import pytest

from onco_dynamics_x.data_engine.survival_metrics import evaluate_ipcw_survival


def test_uncensored_perfect_ranking_has_perfect_time_dependent_auc():
    durations = np.array([1, 2, 3, 4, 5, 6], dtype=float)
    events = np.ones(6, dtype=int)
    predicted = np.array([0.1, 0.2, 0.3, 0.7, 0.8, 0.9])

    report = evaluate_ipcw_survival(
        durations, events, predicted, horizons_months=[3], n_bootstrap=0
    )
    row = report["horizons"][0]

    assert row["metric_status"] == "ok"
    assert row["n_cases"] == 3
    assert row["n_controls"] == 3
    assert row["time_dependent_auc_ipcw"] == pytest.approx(1.0)
    assert row["brier_score_ipcw"] == pytest.approx(0.28 / 6)
    assert len(row["calibration_bins"]) == 6


def test_ipcw_metrics_remain_finite_with_right_censoring():
    durations = np.array([1, 2, 3, 4, 5, 6], dtype=float)
    events = np.array([1, 0, 1, 0, 1, 0], dtype=int)
    predicted = np.array([0.1, 0.25, 0.3, 0.55, 0.7, 0.9])

    report = evaluate_ipcw_survival(
        durations, events, predicted, horizons_months=[3], n_bootstrap=0
    )
    row = report["horizons"][0]

    assert row["metric_status"] == "ok"
    assert row["n_cases"] == 2
    assert row["n_controls"] == 3
    assert 0 <= row["time_dependent_auc_ipcw"] <= 1
    assert row["brier_score_ipcw"] >= 0
    assert 0 <= row["ipcw_observed_survival"] <= 1
    assert row["calibration_gap_observed_minus_predicted"] is not None


def test_horizon_with_no_events_is_reported_as_not_estimable():
    durations = np.array([2, 3, 4, 5], dtype=float)
    events = np.zeros(4, dtype=int)
    predicted = np.array([0.9, 0.8, 0.7, 0.6])

    report = evaluate_ipcw_survival(
        durations, events, predicted, horizons_months=[3], n_bootstrap=0
    )

    assert report["horizons"][0]["metric_status"].startswith("not_estimable")
    assert report["horizons"][0]["time_dependent_auc_ipcw"] is None
    assert report["integrated_brier_score_sparse_grid"] is None


def test_horizons_must_be_sorted_and_unique():
    with pytest.raises(ValueError, match="strictly increasing"):
        evaluate_ipcw_survival(
            [1, 2, 3], [1, 1, 1], [0.1, 0.5, 0.9],
            horizons_months=[36, 12], n_bootstrap=0
        )


def test_bootstrap_is_reproducible():
    durations = np.arange(1, 41, dtype=float)
    events = np.array([1, 0] * 20, dtype=int)
    predicted = np.linspace(0.95, 0.05, 40)

    first = evaluate_ipcw_survival(
        durations, events, predicted,
        horizons_months=[12], n_bootstrap=40, seed=7
    )
    second = evaluate_ipcw_survival(
        durations, events, predicted,
        horizons_months=[12], n_bootstrap=40, seed=7
    )

    assert (
        first["horizons"][0]["bootstrap"]["time_dependent_auc_ipcw_95ci"]
        == second["horizons"][0]["bootstrap"]["time_dependent_auc_ipcw_95ci"]
    )


def test_invalid_prediction_probabilities_are_rejected():
    with pytest.raises(ValueError, match="lie in"):
        evaluate_ipcw_survival(
            [1, 2, 3], [1, 0, 1], [0.2, 1.2, 0.8],
            horizons_months=[2], n_bootstrap=0
        )
