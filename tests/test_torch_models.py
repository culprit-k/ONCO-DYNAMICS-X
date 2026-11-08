import pytest


def _torch_components():
    pytest.importorskip("torch")
    from onco_dynamics_x.ai_engine.multiomics_transformer import ExpressionTransformer
    from onco_dynamics_x.ai_engine.neural_operator import FourierNeuralOperator2D
    from onco_dynamics_x.ai_engine.pinn import FieldPINN, pde_residual, pinn_objective
    from onco_dynamics_x.ai_engine.spatial_gnn import SpatialCellGNN

    return (
        ExpressionTransformer,
        FourierNeuralOperator2D,
        FieldPINN,
        pde_residual,
        pinn_objective,
        SpatialCellGNN,
    )


def test_expression_transformer_output_shape():
    components = _torch_components()
    import torch

    model = components[0](number_of_genes=16, output_features=3, embedding_dim=16)
    assert model(torch.randn(4, 16)).shape == (4, 3)


def test_gnn_returns_one_embedding_per_node():
    components = _torch_components()
    import torch

    model = components[5](input_features=5, hidden_features=12, output_features=2)
    node_features = torch.randn(6, 5)
    edges = torch.tensor([[0, 1, 2, 3, 4], [1, 2, 3, 4, 5]])
    assert model(node_features, edges).shape == (6, 2)


def test_fno_preserves_grid_shape_and_changes_channels():
    components = _torch_components()
    import torch

    model = components[1](
        input_channels=8, output_channels=4, width=12, modes_x=4, modes_y=4, layers=2
    )
    assert model(torch.randn(2, 16, 16, 8)).shape == (2, 16, 16, 4)


def test_pinn_residual_is_differentiable():
    components = _torch_components()
    import torch

    model = components[2](output_fields=2, hidden_width=12, hidden_layers=2)
    coords = torch.rand(10, 3)

    def zero_reaction(values):
        return torch.zeros_like(values)

    residual = components[3](model, coords, [0.1, 0.2], zero_reaction)
    loss, pieces = components[4](model(coords), torch.zeros(10, 2), residual)
    loss.backward()
    assert residual.shape == (10, 2)
    assert torch.isfinite(pieces["total"])
