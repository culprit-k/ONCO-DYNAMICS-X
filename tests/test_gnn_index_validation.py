import pytest


def test_gnn_rejects_negative_edge_indices():
    torch = pytest.importorskip("torch")
    from onco_dynamics_x.ai_engine.spatial_gnn import SpatialCellGNN

    model = SpatialCellGNN(input_features=2, hidden_features=4, output_features=1)
    nodes = torch.randn(3, 2)
    edges = torch.tensor([[0, -1], [1, 2]])
    with pytest.raises(ValueError, match="negative"):
        model(nodes, edges)
