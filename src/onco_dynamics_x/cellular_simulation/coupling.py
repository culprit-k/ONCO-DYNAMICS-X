"""Bidirectional operator-split coupling between PDE fields and cell agents.

The density bridge is a numerical coupling device. Its occupancy scale and
coupling strength are user-specified dimensionless parameters and are not
biologically calibrated.
"""
from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numpy.typing import NDArray

from onco_dynamics_x.cellular_simulation.agent_based import (
    CANCER,
    FIBROBLAST,
    IMMUNE,
    AgentPopulation,
    advance_agents,
)
from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.solver import step
from onco_dynamics_x.state import TumourState

FloatArray = NDArray[np.float64]


def agent_density_fields(
    population: AgentPopulation,
    grid_shape: tuple[int, ...],
    *,
    agents_per_unit: float,
) -> dict[str, FloatArray]:
    """Map discrete agents to per-grid-node density proxies.

    Coordinates use the same nearest-grid convention as the agent/environment
    sampler. agents_per_unit is the number of agents per voxel corresponding
    to a continuum field value of 1.0. It must be selected and documented by
    the experimenter because this module supplies no biological calibration.
    """
    population.validate()
    if not np.isfinite(agents_per_unit) or agents_per_unit <= 0:
        raise ValueError("agents_per_unit must be finite and positive.")
    dimensions = population.spatial_dimensions
    if len(grid_shape) != dimensions or any(size < 2 for size in grid_shape):
        raise ValueError("grid_shape must contain one size >=2 per spatial dimension.")
    counts = {
        "cancer": np.zeros(grid_shape, dtype=float),
        "immune": np.zeros(grid_shape, dtype=float),
        "fibroblast": np.zeros(grid_shape, dtype=float),
    }
    if dimensions == 2:
        row = np.rint((population.y + 1.0) * 0.5 * (grid_shape[0] - 1)).astype(int)
        col = np.rint((population.x + 1.0) * 0.5 * (grid_shape[1] - 1)).astype(int)
        coordinates = (row, col)
    else:
        assert population.z is not None
        ix = np.rint((population.x + 1.0) * 0.5 * (grid_shape[0] - 1)).astype(int)
        iy = np.rint((population.y + 1.0) * 0.5 * (grid_shape[1] - 1)).astype(int)
        iz = np.rint((population.z + 1.0) * 0.5 * (grid_shape[2] - 1)).astype(int)
        coordinates = (ix, iy, iz)

    for phenotype, field_name in (
        (CANCER, "cancer"),
        (IMMUNE, "immune"),
        (FIBROBLAST, "fibroblast"),
    ):
        mask = population.cell_type == phenotype
        np.add.at(counts[field_name], tuple(axis[mask] for axis in coordinates), 1.0)
    return {name: values / agents_per_unit for name, values in counts.items()}


def couple_agents_to_pde(
    continuum_state: TumourState,
    population: AgentPopulation,
    parameters: ModelParameters,
    *,
    coupling_strength: float,
    agents_per_unit: float,
) -> TumourState:
    """Blend PDE cell-density proxies with grid-aggregated agent densities.

    A coupling strength of zero leaves the continuum state unchanged; one
    replaces the three cell-density fields with agent-derived densities.
    Oxygen, nutrient, drug, immune-evasion and vascular fields remain PDE fields.
    """
    continuum_state.validate()
    population.validate()
    parameters.validate()
    if continuum_state.cancer.ndim != population.spatial_dimensions:
        raise ValueError("Agent and PDE environment spatial dimensions do not match.")
    if not np.isfinite(coupling_strength) or not 0.0 <= coupling_strength <= 1.0:
        raise ValueError("coupling_strength must lie in [0, 1].")
    densities = agent_density_fields(
        population, continuum_state.shape, agents_per_unit=agents_per_unit
    )
    result = continuum_state.copy()
    alpha = float(coupling_strength)
    result.cancer = np.clip(
        (1.0 - alpha) * result.cancer + alpha * densities["cancer"],
        0.0,
        parameters.cancer_capacity,
    )
    result.immune = np.clip(
        (1.0 - alpha) * result.immune + alpha * densities["immune"],
        0.0,
        parameters.maximum_field,
    )
    result.fibroblast = np.clip(
        (1.0 - alpha) * result.fibroblast + alpha * densities["fibroblast"],
        0.0,
        parameters.maximum_field,
    )
    result.validate()
    return result


def coupled_step(
    state: TumourState,
    population: AgentPopulation,
    parameters: ModelParameters,
    dx: float,
    dt: float,
    *,
    infusion_rate: float = 0.0,
    coupling_strength: float,
    agents_per_unit: float,
    seed: int = 42,
    maximum_agents: int = 20_000,
) -> tuple[TumourState, AgentPopulation]:
    """Advance PDE fields, agents in the updated microenvironment, then feed agents back."""
    if population.spatial_dimensions != state.cancer.ndim:
        raise ValueError("Agent and PDE environment spatial dimensions must match.")
    pde_updated = step(
        state, parameters, dx=dx, dt=dt, infusion_rate=infusion_rate
    )
    agents_updated = advance_agents(
        population,
        pde_updated,
        parameters,
        dt=dt,
        seed=seed,
        maximum_agents=maximum_agents,
    )
    coupled_state = couple_agents_to_pde(
        pde_updated,
        agents_updated,
        parameters,
        coupling_strength=coupling_strength,
        agents_per_unit=agents_per_unit,
    )
    return coupled_state, agents_updated


def simulate_coupled(
    initial_state: TumourState,
    initial_population: AgentPopulation,
    parameters: ModelParameters,
    dx: float,
    dt: float,
    steps: int,
    *,
    coupling_strength: float,
    agents_per_unit: float,
    infusion_schedule: Callable[[float], float] | None = None,
    save_every: int = 10,
    seed: int = 42,
    maximum_agents: int = 20_000,
) -> tuple[np.ndarray, list[TumourState], list[AgentPopulation]]:
    """Run a reproducible operator-split PDE/ABM simulation with snapshots."""
    if steps < 1 or save_every < 1:
        raise ValueError("steps and save_every must be at least 1.")
    initial_state.validate()
    initial_population.validate()
    if initial_state.cancer.ndim != initial_population.spatial_dimensions:
        raise ValueError("Initial PDE and ABM dimensions must match.")
    if len(initial_population.x) > maximum_agents:
        raise ValueError("Initial population exceeds maximum_agents.")
    if not np.isfinite(coupling_strength) or not 0.0 <= coupling_strength <= 1.0:
        raise ValueError("coupling_strength must lie in [0, 1].")

    current_population = initial_population.copy()
    current_state = couple_agents_to_pde(
        initial_state,
        current_population,
        parameters,
        coupling_strength=coupling_strength,
        agents_per_unit=agents_per_unit,
    )
    state_snapshots = [current_state.copy()]
    population_snapshots = [current_population.copy()]
    times = [0.0]

    for index in range(1, steps + 1):
        time = (index - 1) * dt
        input_rate = (
            0.0 if infusion_schedule is None else float(infusion_schedule(time))
        )
        current_state, current_population = coupled_step(
            current_state,
            current_population,
            parameters,
            dx,
            dt,
            infusion_rate=input_rate,
            coupling_strength=coupling_strength,
            agents_per_unit=agents_per_unit,
            seed=seed + index,
            maximum_agents=maximum_agents,
        )
        if index % save_every == 0 or index == steps:
            times.append(index * dt)
            state_snapshots.append(current_state.copy())
            population_snapshots.append(current_population.copy())

    return (
        np.asarray(times, dtype=float),
        state_snapshots,
        population_snapshots,
    )
