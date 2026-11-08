"""Spatial cell-graph construction from agent coordinates and local PDE fields."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from onco_dynamics_x.cellular_simulation.agent_based import (
    CANCER,
    FIBROBLAST,
    IMMUNE,
    AgentPopulation,
)
from onco_dynamics_x.state import TumourState

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


@dataclass(frozen=True)
class CellGraph:
    """Node feature matrix and directed edge index compatible with message-passing GNNs."""

    node_features: FloatArray
    edge_index: IntArray
    node_ids: tuple[int, ...]

    def validate(self) -> None:
        if self.node_features.ndim != 2:
            raise ValueError("node_features must be two-dimensional.")
        if self.edge_index.ndim != 2 or self.edge_index.shape[0] != 2:
            raise ValueError("edge_index must have shape [2, edges].")
        if self.edge_index.size and (
            int(self.edge_index.min()) < 0
            or int(self.edge_index.max()) >= len(self.node_features)
        ):
            raise ValueError("edge_index contains an out-of-range node index.")
        if len(self.node_ids) != len(self.node_features):
            raise ValueError("node_ids must contain one ID per feature row.")
        if not np.isfinite(self.node_features).all():
            raise ValueError("node_features must be finite.")


def radius_graph(
    population: AgentPopulation,
    radius: float,
    max_neighbors: int | None = 32,
    chunk_size: int = 512,
) -> IntArray:
    """Construct directed edges between cells within a Euclidean radius.

    The chunked calculation avoids allocating a full N-by-N-by-d coordinate
    tensor. With max_neighbors set, each node retains its nearest neighbours.
    """
    population.validate()
    if radius <= 0 or not np.isfinite(radius):
        raise ValueError("radius must be finite and positive.")
    if max_neighbors is not None and max_neighbors < 1:
        raise ValueError("max_neighbors must be positive or None.")
    if chunk_size < 1:
        raise ValueError("chunk_size must be positive.")
    coordinates = (
        np.column_stack([population.x, population.y])
        if population.z is None
        else np.column_stack([population.x, population.y, population.z])
    )
    source_parts: list[IntArray] = []
    target_parts: list[IntArray] = []
    count = len(coordinates)
    radius_squared = radius**2
    all_indices = np.arange(count)

    for start in range(0, count, chunk_size):
        stop = min(start + chunk_size, count)
        differences = coordinates[start:stop, None, :] - coordinates[None, :, :]
        distances_squared = np.sum(differences**2, axis=2)
        for local_index, row in enumerate(distances_squared):
            source_index = start + local_index
            neighbours = np.flatnonzero(
                (row <= radius_squared) & (all_indices != source_index)
            )
            if max_neighbors is not None and len(neighbours) > max_neighbors:
                order = np.argsort(row[neighbours], kind="stable")
                neighbours = neighbours[order[:max_neighbors]]
            if len(neighbours):
                source_parts.append(
                    np.full(len(neighbours), source_index, dtype=np.int64)
                )
                target_parts.append(neighbours.astype(np.int64, copy=False))

    if not source_parts:
        return np.empty((2, 0), dtype=np.int64)
    return np.vstack([np.concatenate(source_parts), np.concatenate(target_parts)])


def _sample_field(field: FloatArray, population: AgentPopulation) -> FloatArray:
    if population.z is None and field.ndim == 2:
        row = np.rint((population.y + 1.0) * 0.5 * (field.shape[0] - 1)).astype(int)
        col = np.rint((population.x + 1.0) * 0.5 * (field.shape[1] - 1)).astype(int)
        return field[row, col]
    if population.z is not None and field.ndim == 3:
        ix = np.rint((population.x + 1.0) * 0.5 * (field.shape[0] - 1)).astype(int)
        iy = np.rint((population.y + 1.0) * 0.5 * (field.shape[1] - 1)).astype(int)
        iz = np.rint((population.z + 1.0) * 0.5 * (field.shape[2] - 1)).astype(int)
        return field[ix, iy, iz]
    raise ValueError("Agent and environment spatial dimensions must match.")


def build_cell_graph(
    population: AgentPopulation,
    environment: TumourState,
    radius: float,
    max_neighbors: int | None = 32,
) -> CellGraph:
    """Build GNN inputs: coordinates, cell-type encoding, age and eight local fields."""
    population.validate()
    environment.validate()
    if population.spatial_dimensions != environment.cancer.ndim:
        raise ValueError("Agent and PDE environment spatial dimensions must match.")
    cell_type_features = np.column_stack(
        [
            population.cell_type == CANCER,
            population.cell_type == IMMUNE,
            population.cell_type == FIBROBLAST,
        ]
    ).astype(float)
    environment_features = np.column_stack(
        [_sample_field(environment.as_dict()[name], population)
         for name in environment.FIELD_NAMES]
    )
    coordinates = (
        np.column_stack([population.x, population.y])
        if population.z is None
        else np.column_stack([population.x, population.y, population.z])
    )
    node_features = np.column_stack(
        [coordinates, cell_type_features, population.age, environment_features]
    ).astype(float)
    graph = CellGraph(
        node_features=node_features,
        edge_index=radius_graph(population, radius, max_neighbors),
        node_ids=tuple(range(len(population.x))),
    )
    graph.validate()
    return graph
