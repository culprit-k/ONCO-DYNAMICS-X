import importlib

import pytest


def _benchmark_modules():
    torch = pytest.importorskip("torch")
    benchmark = importlib.import_module(
        "onco_dynamics_x.ai_engine.synthetic_benchmarks"
    )
    return torch, benchmark


def test_exact_neumann_heat_solution_has_expected_shape_and_nonnegative_values():
    torch, benchmark = _benchmark_modules()
    coordinates = torch.tensor(
        [
            [-1.0, -1.0, 0.0],
            [0.0, 0.0, 0.1],
            [1.0, 1.0, 0.25],
        ],
        dtype=torch.float32,
    )
    values = benchmark.exact_neumann_heat_solution(coordinates, diffusion=0.02)

    assert values.shape == (3, 1)
    assert torch.isfinite(values).all()
    assert torch.all(values >= 0)
    assert values[0, 0].item() == pytest.approx(1.0, abs=1e-6)
    assert values[2, 0] > values[1, 0]


def test_exact_neumann_heat_solution_rejects_negative_time():
    torch, benchmark = _benchmark_modules()
    with pytest.raises(ValueError, match="non-negative"):
        benchmark.exact_neumann_heat_solution(
            torch.tensor([[0.0, 0.0, -0.1]])
        )


def test_periodic_heat_operator_preserves_constant_field():
    torch, benchmark = _benchmark_modules()
    field = torch.full((2, 8, 8), 3.0)
    output = benchmark.periodic_heat_operator(
        field, diffusion=0.05, duration=0.2
    )

    assert output.shape == field.shape
    assert torch.allclose(output, field, atol=1e-6)


def test_periodic_heat_operator_smooths_random_field():
    torch, benchmark = _benchmark_modules()
    generator = torch.Generator().manual_seed(21)
    field = torch.randn((3, 16, 16), generator=generator)
    output = benchmark.periodic_heat_operator(
        field, diffusion=0.05, duration=0.2
    )

    assert output.shape == field.shape
    assert torch.isfinite(output).all()
    assert torch.std(output) < torch.std(field)

def test_synthetic_spatial_graph_is_reproducible_and_has_valid_neighbour_targets():
    torch, benchmark = _benchmark_modules()
    first = benchmark.synthetic_spatial_neighbour_graph(
        num_nodes=48, radius=0.5, seed=7
    )
    second = benchmark.synthetic_spatial_neighbour_graph(
        num_nodes=48, radius=0.5, seed=7
    )
    features, edge_index, targets = first

    assert features.shape == (48, 5)
    assert edge_index.shape[0] == 2
    assert targets.shape == (48, 1)
    assert torch.equal(features, second[0])
    assert torch.equal(edge_index, second[1])
    assert torch.equal(targets, second[2])
    assert torch.all((targets >= 0) & (targets <= 1))
    if edge_index.numel():
        assert torch.all(edge_index[0] != edge_index[1])
