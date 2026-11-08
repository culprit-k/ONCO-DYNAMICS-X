import numpy as np
import pytest

from onco_dynamics_x.calibration.uncertainty import monte_carlo_parameter_uncertainty
from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.solver import diffusion_timestep_limit
from onco_dynamics_x.state import initial_state


def test_parameter_uncertainty_is_reproducible():
    params = ModelParameters()
    state = initial_state(9)
    dx = 2.0 / 8
    dt = min(0.0001, diffusion_timestep_limit(params, dx) * 0.25)
    arguments = {
        "initial": state,
        "parameters": params,
        "dx": dx,
        "dt": dt,
        "steps": 2,
        "parameter_ranges": {"cancer_growth": (0.4, 0.6)},
        "samples": 4,
        "seed": 11,
    }
    first = monte_carlo_parameter_uncertainty(**arguments)
    second = monte_carlo_parameter_uncertainty(**arguments)
    assert np.array_equal(first.final_cancer_mean, second.final_cancer_mean)
    assert first.lower_95 <= first.median <= first.upper_95


def test_parameter_uncertainty_rejects_unknown_parameter():
    with pytest.raises(ValueError, match="Unknown"):
        monte_carlo_parameter_uncertainty(
            initial_state(9), ModelParameters(), 0.25, 0.0001, 1,
            {"not_a_parameter": (0.1, 0.2)}, samples=2
        )
