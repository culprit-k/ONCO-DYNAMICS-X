import numpy as np
import pytest


def _calibration():
    pytest.importorskip("scipy")
    from onco_dynamics_x.calibration.biological import (
        fit_untreated_growth_curve,
        physical_rate_to_dimensionless_rate,
    )

    return fit_untreated_growth_curve, physical_rate_to_dimensionless_rate


def test_logistic_curve_recovers_growth_rate_and_capacity_with_bootstrap():
    fit_untreated_growth_curve, _ = _calibration()
    time = np.arange(0.0, 31.0, 3.0)
    initial, capacity, rate = 120.0, 1_200.0, 0.16
    volume = capacity / (
        1.0 + (capacity / initial - 1.0) * np.exp(-rate * time)
    )

    result = fit_untreated_growth_curve(
        time, volume, n_bootstrap=40, seed=91
    )

    assert result.model == "logistic"
    assert result.parameters["logistic_rate_per_day"] == pytest.approx(
        rate, rel=0.03
    )
    assert result.parameters["carrying_capacity_mm3"] == pytest.approx(
        capacity, rel=0.10
    )
    assert result.locally_identifiable
    assert "logistic_rate_per_day" in result.bootstrap_95ci


def test_exponential_growth_model_is_selected_for_unsaturated_growth():
    fit_untreated_growth_curve, _ = _calibration()
    time = np.arange(0.0, 25.0, 3.0)
    volume = 80.0 * np.exp(0.045 * time)

    result = fit_untreated_growth_curve(time, volume, n_bootstrap=0)

    assert result.model == "exponential"
    assert result.parameters["growth_rate_per_day"] == pytest.approx(
        0.045, rel=1e-5
    )
    assert result.parameters["initial_volume_mm3"] == pytest.approx(80.0)


def test_growth_estimator_rejects_shrinkage_and_nonunique_times():
    fit_untreated_growth_curve, _ = _calibration()
    time = np.arange(5.0)
    shrinkage = 100.0 * np.exp(-0.1 * time)

    with pytest.raises(ValueError, match="positive median growth slope"):
        fit_untreated_growth_curve(time, shrinkage, n_bootstrap=0)
    with pytest.raises(ValueError, match="unique and strictly increasing"):
        fit_untreated_growth_curve(
            [0, 1, 1, 3, 4],
            [100, 110, 120, 140, 160],
            n_bootstrap=0,
        )


def test_growth_estimator_requires_sufficient_longitudinal_observations():
    fit_untreated_growth_curve, _ = _calibration()
    with pytest.raises(ValueError, match="At least five"):
        fit_untreated_growth_curve(
            [0, 1, 2, 3], [100, 110, 125, 140], n_bootstrap=0
        )


def test_physical_to_dimensionless_rate_conversion_is_explicit_and_validated():
    _, convert = _calibration()
    assert convert(0.12, 7.0) == pytest.approx(0.84)
    with pytest.raises(ValueError, match="positive"):
        convert(-0.1, 2.0)
