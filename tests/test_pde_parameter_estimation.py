from dataclasses import replace

import numpy as np
import pytest


def _estimator():
    pytest.importorskip("scipy")
    from onco_dynamics_x.calibration.pde_parameter_estimation import (
        fit_pde_parameter_subset,
    )

    return fit_pde_parameter_subset


def _synthetic_observations(parameters, initial, times, dt, dx, fields):
    from onco_dynamics_x.solver import step

    current = initial.copy()
    result = {name: [getattr(current, name).copy()] for name in fields}
    for index in range(1, len(times)):
        interval_steps = int(round((times[index] - times[index - 1]) / dt))
        for _ in range(interval_steps):
            current = step(current, parameters, dx, dt)
        for name in fields:
            result[name].append(getattr(current, name).copy())
    return {name: np.stack(values) for name, values in result.items()}


def test_pde_parameter_estimator_recovers_one_known_synthetic_coefficient():
    from onco_dynamics_x.parameters import ModelParameters
    from onco_dynamics_x.state import initial_state

    fit_pde_parameter_subset = _estimator()
    initial = initial_state(grid_size=7)
    baseline = ModelParameters()
    truth = replace(baseline, cancer_growth=0.72)
    dt, dx = 0.01, 0.25
    times = np.array([0.0, 0.1, 0.2, 0.3])
    observations = _synthetic_observations(
        truth, initial, times, dt, dx, ("cancer", "oxygen")
    )

    result = fit_pde_parameter_subset(
        times,
        observations,
        initial,
        baseline,
        parameter_names=("cancer_growth",),
        parameter_bounds={"cancer_growth": (0.2, 1.2)},
        observation_sd={"cancer": 0.01, "oxygen": 0.01},
        dx=dx,
        dt=dt,
    )

    assert result.optimiser_success
    assert result.locally_identifiable
    assert result.parameter_estimates["cancer_growth"] == pytest.approx(
        truth.cancer_growth, rel=0.03
    )
    assert result.parameter_standard_errors is not None
    assert result.observation_count == 4 * (7 * 7) * 2


def test_pde_parameter_estimator_rejects_bad_observation_alignment_and_bounds():
    from onco_dynamics_x.parameters import ModelParameters
    from onco_dynamics_x.state import initial_state

    fit_pde_parameter_subset = _estimator()
    initial = initial_state(grid_size=5)
    times = [0.0, 0.1]
    observations = {
        "cancer": np.stack([initial.cancer, initial.cancer]),
    }
    with pytest.raises(ValueError, match="time zero must match"):
        altered = {"cancer": observations["cancer"].copy()}
        altered["cancer"][0, 0, 0] += 0.1
        fit_pde_parameter_subset(
            times,
            altered,
            initial,
            ModelParameters(),
            parameter_names=("cancer_growth",),
            parameter_bounds={"cancer_growth": (0.1, 1.0)},
            observation_sd={"cancer": 0.1},
            dx=0.5,
            dt=0.01,
        )
    with pytest.raises(ValueError, match="positive and ordered"):
        fit_pde_parameter_subset(
            times,
            observations,
            initial,
            ModelParameters(),
            parameter_names=("cancer_growth",),
            parameter_bounds={"cancer_growth": (0.0, 1.0)},
            observation_sd={"cancer": 0.1},
            dx=0.5,
            dt=0.01,
        )


def test_pde_parameter_estimator_limits_each_stage_to_three_parameters():
    from onco_dynamics_x.parameters import ModelParameters
    from onco_dynamics_x.state import initial_state

    fit_pde_parameter_subset = _estimator()
    initial = initial_state(grid_size=5)
    observations = {"cancer": np.stack([initial.cancer, initial.cancer])}
    with pytest.raises(ValueError, match="one to three"):
        fit_pde_parameter_subset(
            [0.0, 0.1],
            observations,
            initial,
            ModelParameters(),
            parameter_names=(
                "cancer_growth",
                "baseline_cancer_loss",
                "immune_killing",
                "oxygen_consumption",
            ),
            parameter_bounds={
                "cancer_growth": (0.1, 1.0),
                "baseline_cancer_loss": (0.001, 0.1),
                "immune_killing": (0.1, 1.0),
                "oxygen_consumption": (0.01, 1.0),
            },
            observation_sd={"cancer": 0.1},
            dx=0.5,
            dt=0.01,
        )

def test_pde_parameter_estimator_withholds_covariance_for_unobservable_drug_decay():
    from onco_dynamics_x.parameters import ModelParameters
    from onco_dynamics_x.state import initial_state

    fit_pde_parameter_subset = _estimator()
    initial = initial_state(grid_size=5)
    assert np.all(initial.drug == 0.0)
    times = [0.0, 0.1, 0.2]
    # A zero-drug experiment contains no information about drug-field decay.
    observations = {
        "cancer": np.stack([initial.cancer, initial.cancer, initial.cancer]),
    }

    result = fit_pde_parameter_subset(
        times,
        observations,
        initial,
        ModelParameters(),
        parameter_names=("drug_decay",),
        parameter_bounds={"drug_decay": (0.01, 1.0)},
        observation_sd={"cancer": 0.05},
        dx=0.5,
        dt=0.01,
    )

    assert not result.locally_identifiable
    assert result.jacobian_rank == 0
    assert result.parameter_covariance is None
    assert result.parameter_standard_errors is None
    assert result.warnings
