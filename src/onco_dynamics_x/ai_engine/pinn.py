"""Optional PyTorch physics-informed neural network utilities."""
from __future__ import annotations

from collections.abc import Callable, Sequence

import torch
from torch import Tensor, nn


class FieldPINN(nn.Module):
    """MLP mapping (x, y, t) coordinates to eight non-negative model fields."""

    def __init__(
        self,
        output_fields: int = 8,
        hidden_width: int = 128,
        hidden_layers: int = 5,
    ) -> None:
        super().__init__()
        if output_fields < 1 or hidden_width < 4 or hidden_layers < 1:
            raise ValueError("Network dimensions must be positive and non-trivial.")
        modules: list[nn.Module] = [nn.Linear(3, hidden_width), nn.Tanh()]
        for _ in range(hidden_layers - 1):
            modules.extend([nn.Linear(hidden_width, hidden_width), nn.Tanh()])
        modules.append(nn.Linear(hidden_width, output_fields))
        self.network = nn.Sequential(*modules)
        self.positive_output = nn.Softplus()

    def forward(self, coordinates: Tensor) -> Tensor:
        if coordinates.ndim != 2 or coordinates.shape[1] != 3:
            raise ValueError("coordinates must have shape [batch, 3] for (x, y, t).")
        return self.positive_output(self.network(coordinates))


def pde_residual(
    model: nn.Module,
    coordinates: Tensor,
    diffusion_coefficients: Sequence[float],
    reaction_function: Callable[[Tensor], Tensor],
) -> Tensor:
    """Compute u_t - D*(u_xx+u_yy) - R(u) using automatic differentiation.

    The reaction function receives the predicted fields [batch, n_fields] and
    returns a tensor of the same shape. Coordinates are ordered (x, y, t).
    """
    if not coordinates.requires_grad:
        coordinates = coordinates.clone().detach().requires_grad_(True)
    prediction = model(coordinates)
    if prediction.ndim != 2 or prediction.shape[1] != len(diffusion_coefficients):
        raise ValueError("Diffusion coefficients must match the model output fields.")
    reactions = reaction_function(prediction)
    if reactions.shape != prediction.shape:
        raise ValueError("reaction_function must return the same shape as model output.")

    residuals = []
    for field_index, diffusion in enumerate(diffusion_coefficients):
        if diffusion < 0:
            raise ValueError("Diffusion coefficients must be non-negative.")
        first = torch.autograd.grad(
            prediction[:, field_index].sum(),
            coordinates,
            create_graph=True,
            retain_graph=True,
        )[0]
        second_x = torch.autograd.grad(
            first[:, 0].sum(), coordinates, create_graph=True, retain_graph=True
        )[0][:, 0]
        second_y = torch.autograd.grad(
            first[:, 1].sum(), coordinates, create_graph=True, retain_graph=True
        )[0][:, 1]
        residual = (
            first[:, 2]
            - float(diffusion) * (second_x + second_y)
            - reactions[:, field_index]
        )
        residuals.append(residual)
    return torch.stack(residuals, dim=1)


def pinn_objective(
    data_predictions: Tensor,
    data_targets: Tensor,
    physics_residuals: Tensor,
    *,
    data_weight: float = 1.0,
    physics_weight: float = 1.0,
) -> tuple[Tensor, dict[str, Tensor]]:
    """Combine observation and PDE residual losses with explicit weights."""
    if data_predictions.shape != data_targets.shape:
        raise ValueError("Data predictions and targets must have matching shapes.")
    if data_weight < 0 or physics_weight < 0:
        raise ValueError("Loss weights must be non-negative.")
    data_loss = torch.mean((data_predictions - data_targets) ** 2)
    physics_loss = torch.mean(physics_residuals**2)
    total = data_weight * data_loss + physics_weight * physics_loss
    return total, {"data": data_loss, "physics": physics_loss, "total": total}
