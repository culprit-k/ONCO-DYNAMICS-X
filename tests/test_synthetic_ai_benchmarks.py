import pytest

torch = pytest.importorskip("torch")

from onco_dynamics_x.ai_engine.synthetic_benchmarks import (
    exact_neumann_heat_solution,
    periodic_heat_operator,
)


def test_exact_neumann_heat_solution_has_expected_shape_and_nonnegative_values():
    coordinates = torch.tensor(
        [
            [-1.0, -1.0, 0.0],
            [0.0, 0.0, 0.1],
            [1.0, 1.0, 0.25],
        ],
        dtype=torch.float32,
    )
    values = exact_neumann_heat_solution(coordinates, diffusion=0.02)

    assert values.shape == (3, 1)
    assert torch.isfinite(values).all()
    assert torch.all(values >= 0)
    assert values[0, 0].item() == pytest.approx(0.0, abs=1e-6)
    assert values[1, 0] > values[2, 0]


def test_exact_neumann_heat_solution_rejects_negative_time():
    with pytest.raises(ValueError, match="non-negative"):
        exact_neumann_heat_solution(torch.tensor([[0.0, 0.0, -0.1]]))


def test_periodic_heat_operator_preserves_constant_field():
    field = torch.full((2, 8, 8), 3.0)
    output = periodic_heat_operator(field, diffusion=0.05, duration=0.2)

    assert output.shape == field.shape
    assert torch.allclose(output, field, atol=1e-6)


def test_periodic_heat_operator_smooths_random_field():
    generator = torch.Generator().manual_seed(21)
    field = torch.randn((3, 16, 16), generator=generator)
    output = periodic_heat_operator(field, diffusion=0.05, duration=0.2)

    assert output.shape == field.shape
    assert torch.isfinite(output).all()
    assert torch.std(output) < torch.std(field)
