import pytest

from onco_dynamics_x.cellular_simulation.agent_based import advance_agents, initialize_agents
from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.state import initial_state


def test_agent_model_rejects_3d_environment_until_3d_agent_rules_exist():
    population = initialize_agents(cancer_count=2, immune_count=2, fibroblast_count=1)
    environment = initial_state(grid_size=7, spatial_dimensions=3)
    with pytest.raises(ValueError, match="requires a 2D environment"):
        advance_agents(population, environment, ModelParameters(), dt=0.001)
