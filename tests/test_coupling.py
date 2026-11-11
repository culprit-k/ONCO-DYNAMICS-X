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


def test_simulate_coupled_returns_reproducible_snapshots_and_time_grid():
    from onco_dynamics_x.cellular_simulation.coupling import simulate_coupled

    parameters = ModelParameters()
    state = initial_state(7)
    population = initialize_agents(6, 4, 3, seed=18)
    dx = 2.0 / 6
    dt = min(1e-5, diffusion_timestep_limit(parameters, dx) * 0.1)
    kwargs = {
        "coupling_strength": 0.2,
        "agents_per_unit": 2.0,
        "save_every": 2,
        "seed": 91,
    }

    first = simulate_coupled(
        state, population, parameters, dx, dt, steps=3, **kwargs
    )
    second = simulate_coupled(
        state, population, parameters, dx, dt, steps=3, **kwargs
    )
    first_times, first_states, first_populations = first
    second_times, second_states, second_populations = second

    assert np.allclose(first_times, [0.0, 2.0 * dt, 3.0 * dt])
    assert len(first_states) == 3
    assert len(first_populations) == 3
    assert np.array_equal(first_times, second_times)
    for left, right in zip(first_states, second_states, strict=True):
        left.validate()
        for name in left.FIELD_NAMES:
            assert np.allclose(left.as_dict()[name], right.as_dict()[name])
    for left, right in zip(first_populations, second_populations, strict=True):
        left.validate()
        assert np.array_equal(left.cell_type, right.cell_type)
        assert np.array_equal(left.x, right.x)


def test_cloud_in_cell_deposition_preserves_agent_counts_and_spreads_mass_in_2d():
    population = AgentPopulation(
        x=np.array([0.25, 0.25, 0.25]),
        y=np.array([0.25, 0.25, 0.25]),
        cell_type=np.array([CANCER, IMMUNE, FIBROBLAST], dtype=np.int64),
        age=np.zeros(3, dtype=float),
    )
    fields = agent_density_fields(
        population,
        (5, 5),
        agents_per_unit=1.0,
        deposition_method="cloud_in_cell",
    )

    for name in ("cancer", "immune", "fibroblast"):
        assert np.sum(fields[name]) == pytest.approx(1.0)
        assert np.count_nonzero(fields[name]) == 4
        assert np.allclose(fields[name][2:4, 2:4], 0.25)


def test_cloud_in_cell_deposition_preserves_agent_counts_in_3d():
    population = AgentPopulation(
        x=np.array([0.25, -0.25, 0.25]),
        y=np.array([0.25, 0.25, -0.25]),
        z=np.array([0.25, 0.25, 0.25]),
        cell_type=np.array([CANCER, IMMUNE, FIBROBLAST], dtype=np.int64),
        age=np.zeros(3, dtype=float),
    )
    fields = agent_density_fields(
        population,
        (5, 5, 5),
        agents_per_unit=2.0,
        deposition_method="cloud_in_cell",
    )

    for name in ("cancer", "immune", "fibroblast"):
        assert np.sum(fields[name]) == pytest.approx(0.5)
        assert np.count_nonzero(fields[name]) == 8


def test_deposition_methods_are_validated_and_nearest_grid_remains_default():
    population = _three_cell_population_2d()
    nearest = agent_density_fields(
        population, (5, 5), agents_per_unit=1.0
    )
    explicit_nearest = agent_density_fields(
        population,
        (5, 5),
        agents_per_unit=1.0,
        deposition_method="nearest_grid",
    )
    assert np.array_equal(nearest["cancer"], explicit_nearest["cancer"])
    with pytest.raises(ValueError, match="deposition_method"):
        agent_density_fields(
            population,
            (5, 5),
            agents_per_unit=1.0,
            deposition_method="kernel",
        )


def test_coupled_step_accepts_cloud_in_cell_deposition():
    state = initial_state(7)
    population = initialize_agents(6, 4, 3, seed=18)
    parameters = ModelParameters()
    dx = 2.0 / 6
    dt = min(1e-5, diffusion_timestep_limit(parameters, dx) * 0.1)
    result_state, result_population = coupled_step(
        state,
        population,
        parameters,
        dx,
        dt,
        coupling_strength=0.2,
        agents_per_unit=2.0,
        deposition_method="cloud_in_cell",
        seed=91,
    )

    result_state.validate()
    result_population.validate()
