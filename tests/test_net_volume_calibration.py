import numpy as np
import pytest


def _net_rate():
    pytest.importorskip("scipy")
    from onco_dynamics_x.calibration.net_volume import (
        bootstrap_group_net_rates,
        fit_treatment_exposed_net_volume_rate,
    )

    return fit_treatment_exposed_net_volume_rate, bootstrap_group_net_rates


def test_net_volume_rate_recovers_treatment_exposed_shrinkage_slope():
    fit_rate, _ = _net_rate()
    days = np.array([0.0, 21.0, 70.0, 140.0])
    volume = 900.0 * np.exp(-0.012 * days)

    fit = fit_rate(days, volume)

    assert fit.net_log_rate_per_day == pytest.approx(-0.012, rel=1e-10)
    assert fit.halving_time_days == pytest.approx(np.log(2.0) / 0.012)
    assert fit.doubling_time_days is None
    assert fit.r_squared == pytest.approx(1.0)
    assert fit.n_observations == 4


def test_net_volume_rate_recovers_growth_and_does_not_assume_shrinkage():
    fit_rate, _ = _net_rate()
    days = [0.0, 10.0, 30.0]
    volume = 150.0 * np.exp(0.02 * np.asarray(days))

    fit = fit_rate(days, volume)

    assert fit.net_log_rate_per_day == pytest.approx(0.02, rel=1e-10)
    assert fit.doubling_time_days == pytest.approx(np.log(2.0) / 0.02)
    assert fit.halving_time_days is None


def test_net_volume_rate_rejects_censored_zero_values_and_duplicate_days():
    fit_rate, _ = _net_rate()
    with pytest.raises(ValueError, match="left-censored"):
        fit_rate([0, 21, 42], [100, 50, 0])
    with pytest.raises(ValueError, match="unique and strictly increasing"):
        fit_rate([0, 21, 21], [100, 80, 60])


def test_net_volume_rate_requires_three_visits_for_uncertainty():
    fit_rate, _ = _net_rate()
    with pytest.raises(ValueError, match="At least three"):
        fit_rate([0, 21], [100, 80])


def test_group_bootstrap_uses_subject_level_rates():
    _, bootstrap = _net_rate()
    result = bootstrap([-0.01, -0.02, 0.005, -0.015], n_bootstrap=500, seed=7)

    assert result["subject_series_count"] == 4
    assert result["median_net_log_rate_per_day"] == pytest.approx(-0.0125)
    assert len(result["median_rate_bootstrap_95ci_per_day"]) == 2
