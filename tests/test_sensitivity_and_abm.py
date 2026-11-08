import numpy as np

from onco_dynamics_x.calibration.sensitivity import local_parameter_sensitivity
from onco_dynamics_x.cellular_simulation.agent_based import (
    CANCER,
    IMMUNE,
    advance_agents,
    initialize_agents,
)
from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.solver import diffusion_timestep_limit
from onco_dynamics_x.state import initial_state


def test_agent_initialization_is_reproducible():
    first = initialize_agents(seed=9)
    second = initialize_agents(seed=9)
    assert np.array_equal(first.x, second.x)
    assert np.array_equal(first.cell_type, second.cell_type)


def test_agent_model_preserves_coordinates_and_known_phenotypes():
    environment = initial_state(21)
    population = initialize_agents(20, 12, 8, seed=3)
    result = advance_agents(population, environment, ModelParameters(), dt=0.02, seed=5)
    result.validate()
    assert result.count(CANCER) >= 0
    assert result.count(IMMUNE) >= 0
    assert np.all((result.x >= -1) & (result.x <= 1))
    assert np.all((result.y >= -1) & (result.y <= 1))


def test_sensitivity_returns_finite_elasticity():
    params = ModelParameters()
    state = initial_state(9)
    dx = 2.0 / 8
    dt = min(0.0001, diffusion_timestep_limit(params, dx) * 0.25)
    result = local_parameter_sensitivity(
        "cancer_growth", params, state, dx, dt, steps=2
    )
    assert np.isfinite(result["elasticity"])
    assert result["parameter_value"] == params.cancer_growth
