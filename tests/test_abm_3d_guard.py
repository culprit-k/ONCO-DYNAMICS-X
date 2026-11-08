import numpy as np

from onco_dynamics_x.cellular_simulation.agent_based import advance_agents, initialize_agents
from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.state import initial_state


def test_agent_model_runs_in_3d_environment():
    population = initialize_agents(
        cancer_count=3,
        immune_count=2,
        fibroblast_count=1,
        seed=4,
        spatial_dimensions=3,
    )
    environment = initial_state(grid_size=7, spatial_dimensions=3)
    updated = advance_agents(
        population, environment, ModelParameters(), dt=0.001, seed=10
    )
    updated.validate()
    assert updated.spatial_dimensions == 3
    assert updated.z is not None
    assert np.all((-1.0 <= updated.z) & (updated.z <= 1.0))
