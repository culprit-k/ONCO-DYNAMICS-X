import numpy as np
import pytest

from onco_dynamics_x.cellular_simulation.agent_based import (
    CANCER,
    IMMUNE,
    advance_agents,
    initialize_agents,
)
from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.state import initial_state


def test_3d_agent_initialization_is_reproducible():
    first = initialize_agents(20, 12, 8, seed=4, spatial_dimensions=3)
    second = initialize_agents(20, 12, 8, seed=4, spatial_dimensions=3)
    assert first.z is not None
    assert np.array_equal(first.z, second.z)
    assert first.spatial_dimensions == 3


def test_3d_agent_model_preserves_coordinates_and_cell_types():
    environment = initial_state(grid_size=9, spatial_dimensions=3)
    population = initialize_agents(20, 12, 8, seed=3, spatial_dimensions=3)
    result = advance_agents(
        population, environment, ModelParameters(), dt=0.02, seed=5
    )
    result.validate()
    assert result.z is not None
    assert np.all((-1 <= result.z) & (result.z <= 1))
    assert result.count(CANCER) >= 0
    assert result.count(IMMUNE) >= 0


def test_agent_and_environment_dimensions_must_match():
    population = initialize_agents(2, 2, 1, spatial_dimensions=2)
    environment = initial_state(grid_size=7, spatial_dimensions=3)
    with pytest.raises(ValueError, match="dimensions must match"):
        advance_agents(population, environment, ModelParameters(), dt=0.01)
