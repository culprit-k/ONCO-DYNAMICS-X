import numpy as np
import pytest

from onco_dynamics_x.cellular_simulation.agent_based import (
    CANCER,
    FIBROBLAST,
    IMMUNE,
    AgentPopulation,
    initialize_agents,
)
from onco_dynamics_x.cellular_simulation.coupling import (
    agent_density_fields,
    couple_agents_to_pde,
    coupled_step,
)
from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.solver import diffusion_timestep_limit
from onco_dynamics_x.state import initial_state


def _three_cell_population_2d():
    return AgentPopulation(
        x=np.array([0.0, 0.0, 0.0]),
        y=np.array([0.0, 0.0, 0.0]),
        cell_type=np.array([CANCER, IMMUNE, FIBROBLAST], dtype=np.int64),
        age=np.zeros(3, dtype=float),
    )


def test_agent_density_fields_accumulate_each_phenotype_in_2d():
    population = _three_cell_population_2d()
    fields = agent_density_fields(population, (5, 5), agents_per_unit=2.0)

    assert fields["cancer"][2, 2] == pytest.approx(0.5)
    assert fields["immune"][2, 2] == pytest.approx(0.5)
    assert fields["fibroblast"][2, 2] == pytest.approx(0.5)
    assert np.sum(fields["cancer"]) == pytest.approx(0.5)


def test_agent_density_fields_support_3d_coordinates():
    population = AgentPopulation(
        x=np.array([0.0]),
        y=np.array([0.0]),
        z=np.array([0.0]),
        cell_type=np.array([CANCER], dtype=np.int64),
        age=np.zeros(1, dtype=float),
    )
    fields = agent_density_fields(population, (5, 5, 5), agents_per_unit=1.0)

    assert fields["cancer"][2, 2, 2] == pytest.approx(1.0)
    assert np.count_nonzero(fields["cancer"]) == 1


def test_zero_coupling_preserves_continuum_fields():
    state = initial_state(9)
    population = _three_cell_population_2d()
    result = couple_agents_to_pde(
        state,
        population,
        ModelParameters(),
        coupling_strength=0.0,
        agents_per_unit=1.0,
    )

    for name in state.FIELD_NAMES:
        assert np.array_equal(result.as_dict()[name], state.as_dict()[name])


def test_full_coupling_replaces_cell_fields_but_not_microenvironment_fields():
    state = initial_state(5)
    population = _three_cell_population_2d()
    result = couple_agents_to_pde(
        state,
        population,
        ModelParameters(),
        coupling_strength=1.0,
        agents_per_unit=2.0,
    )

    assert result.cancer[2, 2] == pytest.approx(0.5)
    assert result.immune[2, 2] == pytest.approx(0.5)
    assert result.fibroblast[2, 2] == pytest.approx(0.5)
    for name in ("oxygen", "nutrient", "drug", "evasion", "vasculature"):
        assert np.array_equal(result.as_dict()[name], state.as_dict()[name])


def test_coupled_step_is_reproducible_and_preserves_valid_state():
    state = initial_state(11)
    population = initialize_agents(8, 5, 4, seed=12)
    parameters = ModelParameters()
    dx = 2.0 / 10
    dt = min(1e-5, diffusion_timestep_limit(parameters, dx) * 0.25)
    kwargs = {
        "coupling_strength": 0.15,
        "agents_per_unit": 2.0,
        "seed": 77,
    }

    first_state, first_population = coupled_step(
        state, population, parameters, dx, dt, **kwargs
    )
    second_state, second_population = coupled_step(
        state, population, parameters, dx, dt, **kwargs
    )

    first_state.validate()
    first_population.validate()
    for name in first_state.FIELD_NAMES:
        assert np.allclose(
            first_state.as_dict()[name], second_state.as_dict()[name]
        )
    assert np.array_equal(first_population.x, second_population.x)
    assert np.array_equal(first_population.cell_type, second_population.cell_type)


def test_invalid_density_scaling_and_coupling_are_rejected():
    population = _three_cell_population_2d()
    with pytest.raises(ValueError, match="agents_per_unit"):
        agent_density_fields(population, (5, 5), agents_per_unit=0.0)
    with pytest.raises(ValueError, match="coupling_strength"):
        couple_agents_to_pde(
            initial_state(5),
            population,
            ModelParameters(),
            coupling_strength=1.5,
            agents_per_unit=1.0,
        )
