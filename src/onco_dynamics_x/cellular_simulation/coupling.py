"""Bidirectional operator-split coupling between PDE fields and cell agents.

The density bridge is a numerical coupling device. Its occupancy scale and
coupling strength are user-specified dimensionless parameters and are not
biologically calibrated.
"""
from __future__ import annotations

from collections.abc import Callable
from itertools import product
from typing import Literal

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
DepositionMethod = Literal["nearest_grid", "cloud_in_cell"]


def _continuous_grid_coordinates(
    population: AgentPopulation,
    grid_shape: tuple[int, ...],
) -> tuple[FloatArray, ...]:
    """Convert model coordinates to continuous grid indices in axis order."""
    if population.spatial_dimensions == 2:
        return (
            (population.y + 1.0) * 0.5 * (grid_shape[0] - 1),
            (population.x + 1.0) * 0.5 * (grid_shape[1] - 1),
        )
    assert population.z is not None
    return (
        (population.x + 1.0) * 0.5 * (grid_shape[0] - 1),
        (population.y + 1.0) * 0.5 * (grid_shape[1] - 1),
        (population.z + 1.0) * 0.5 * (grid_shape[2] - 1),
    )


def _deposit_nearest(
    field: FloatArray,
    coordinates: tuple[FloatArray, ...],
    mask: NDArray[np.bool_],
) -> None:
    indices = tuple(np.rint(axis[mask]).astype(int) for axis in coordinates)
    np.add.at(field, indices, 1.0)


def _deposit_cloud_in_cell(
    field: FloatArray,
    coordinates: tuple[FloatArray, ...],
    mask: NDArray[np.bool_],
) -> None:
    """Deposit agents to adjacent grid nodes with multilinear weights.

    For in-domain coordinates, weights for each agent sum to one, preserving
    the total deposited agent count before any later continuum field clipping.
    """
    if not np.any(mask):
        return
    lower = [np.floor(axis[mask]).astype(int) for axis in coordinates]
    fractions = [
        axis[mask] - low
        for axis, low in zip(coordinates, lower, strict=True)
    ]
    upper_limits = np.asarray(field.shape, dtype=int) - 1

    for corner in product((0, 1), repeat=field.ndim):
        indices = []
        weights = np.ones(int(np.count_nonzero(mask)), dtype=float)
        for axis_index, bit in enumerate(corner):
            if bit:
                indices.append(
                    np.minimum(lower[axis_index] + 1, upper_limits[axis_index])
                )
                weights *= fractions[axis_index]
            else:
                indices.append(lower[axis_index])
                weights *= 1.0 - fractions[axis_index]
        np.add.at(field, tuple(indices), weights)


def agent_density_fields(
    population: AgentPopulation,
    grid_shape: tuple[int, ...],
    *,
    agents_per_unit: float,
    deposition_method: DepositionMethod = "nearest_grid",
) -> dict[str, FloatArray]:
    """Map discrete agents to per-grid-node density proxies.

    nearest_grid assigns each agent to one node. cloud_in_cell uses multilinear
    weights across adjacent nodes and conserves the total deposited count before
    any continuum field clipping. Neither method supplies biological calibration;
    agents_per_unit remains an explicit dimensionless occupancy scale.
    """
    population.validate()
    if not np.isfinite(agents_per_unit) or agents_per_unit <= 0:
        raise ValueError("agents_per_unit must be finite and positive.")
    if deposition_method not in {"nearest_grid", "cloud_in_cell"}:
        raise ValueError(
            "deposition_method must be 'nearest_grid' or 'cloud_in_cell'."
        )
    dimensions = population.spatial_dimensions
    if len(grid_shape) != dimensions or any(size < 2 for size in grid_shape):
        raise ValueError("grid_shape must contain one size >=2 per spatial dimension.")

    counts = {
        "cancer": np.zeros(grid_shape, dtype=float),
        "immune": np.zeros(grid_shape, dtype=float),
        "fibroblast": np.zeros(grid_shape, dtype=float),
    }
    coordinates = _continuous_grid_coordinates(population, grid_shape)
    for phenotype, field_name in (
        (CANCER, "cancer"),
        (IMMUNE, "immune"),
        (FIBROBLAST, "fibroblast"),
    ):
        mask = population.cell_type == phenotype
        if deposition_method == "nearest_grid":
            _deposit_nearest(counts[field_name], coordinates, mask)
        else:
            _deposit_cloud_in_cell(counts[field_name], coordinates, mask)
    return {name: values / agents_per_unit for name, values in counts.items()}


def couple_agents_to_pde(
    continuum_state: TumourState,
    population: AgentPopulation,
    parameters: ModelParameters,
    *,
    coupling_strength: float,
    agents_per_unit: float,
    deposition_method: DepositionMethod = "nearest_grid",
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
        population,
        continuum_state.shape,
        agents_per_unit=agents_per_unit,
        deposition_method=deposition_method,
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
    deposition_method: DepositionMethod = "nearest_grid",
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
        deposition_method=deposition_method,
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
    deposition_method: DepositionMethod = "nearest_grid",
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
    if deposition_method not in {"nearest_grid", "cloud_in_cell"}:
        raise ValueError(
            "deposition_method must be 'nearest_grid' or 'cloud_in_cell'."
        )

    current_population = initial_population.copy()
    current_state = couple_agents_to_pde(
        initial_state,
        current_population,
        parameters,
        coupling_strength=coupling_strength,
        agents_per_unit=agents_per_unit,
        deposition_method=deposition_method,
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
            deposition_method=deposition_method,
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
