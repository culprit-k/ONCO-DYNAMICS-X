import numpy as np
import pytest

from onco_dynamics_x.cellular_simulation.agent_based import initialize_agents
from onco_dynamics_x.cellular_simulation.graph import build_cell_graph, radius_graph
from onco_dynamics_x.data_engine.splits import group_train_test_split
from onco_dynamics_x.state import initial_state


def test_radius_graph_creates_nearby_directed_edges_without_self_loops():
    population = initialize_agents(4, 0, 0, seed=1)
    population.x[:] = np.array([0.0, 0.05, 0.8, -0.7])
    population.y[:] = 0.0
    edges = radius_graph(population, radius=0.1, max_neighbors=None)
    assert edges.shape[0] == 2
    assert not np.any(edges[0] == edges[1])
    assert (0, 1) in set(map(tuple, edges.T.tolist()))
    assert (1, 0) in set(map(tuple, edges.T.tolist()))
    assert not any(2 in pair for pair in edges.T.tolist())


def test_cell_graph_has_expected_features_for_2d_and_3d():
    for dimensions in (2, 3):
        population = initialize_agents(5, 3, 2, seed=6, spatial_dimensions=dimensions)
        environment = initial_state(grid_size=7, spatial_dimensions=dimensions)
        graph = build_cell_graph(population, environment, radius=0.4)
        graph.validate()
        expected_features = dimensions + 3 + 1 + 8
        assert graph.node_features.shape == (10, expected_features)
        assert graph.edge_index.shape[0] == 2


def test_radius_graph_rejects_invalid_neighbor_count():
    population = initialize_agents(2, 0, 0)
    with pytest.raises(ValueError, match="max_neighbors"):
        radius_graph(population, radius=0.5, max_neighbors=0)


def test_group_split_keeps_each_patient_in_only_one_partition():
    patient_ids = ["P1", "P1", "P2", "P3", "P3", "P4", "P5", "P5"]
    train, test = group_train_test_split(patient_ids, test_size=0.4, random_seed=5)
    train_patients = {patient_ids[index] for index in train}
    test_patients = {patient_ids[index] for index in test}
    assert train_patients.isdisjoint(test_patients)
    assert len(train) + len(test) == len(patient_ids)


def test_group_split_is_reproducible_and_rejects_missing_groups():
    groups = ["P1", "P2", "P2", "P3", "P4"]
    first = group_train_test_split(groups, random_seed=12)
    second = group_train_test_split(groups, random_seed=12)
    assert np.array_equal(first[0], second[0])
    assert np.array_equal(first[1], second[1])
    with pytest.raises(ValueError, match="empty"):
        group_train_test_split(["P1", "", "P2"])
    with pytest.raises(ValueError, match="empty"):
        group_train_test_split(["P1", np.nan, "P2"])
