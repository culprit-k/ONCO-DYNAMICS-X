"""Analytic synthetic operators used to benchmark optional AI architectures.

These operators provide known mathematical ground truth for software tests.
They are not breast-cancer models or patient-derived training data.
"""
from __future__ import annotations

import math

import torch
from torch import Tensor


def exact_neumann_heat_solution(
    coordinates: Tensor,
    diffusion: float = 0.02,
) -> Tensor:
    """Exact non-negative heat-equation solution on [-1,1]^2 with no-flux edges.

    u(x,y,t) = h(x,t) * h(y,t), where
    h(q,t) = (1 + exp(-D*pi^2*t)*cos(pi*q))/2.
    It is non-negative, has zero normal derivative at x/y = +/-1,
    and satisfies u_t - D*(u_xx+u_yy)=0.
    Coordinates are [batch,3] ordered x,y,t.
    """
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("coordinates must have shape [batch, 3] for x, y, t.")
    if diffusion <= 0 or not math.isfinite(diffusion):
        raise ValueError("diffusion must be finite and positive.")
    x, y, time = coordinates.unbind(dim=1)
    if torch.any(time < 0):
        raise ValueError("time coordinates must be non-negative.")
    decay = torch.exp(-diffusion * math.pi**2 * time)
    field_x = 0.5 * (1.0 - decay * torch.cos(math.pi * x))
    field_y = 0.5 * (1.0 - decay * torch.cos(math.pi * y))
    return (field_x * field_y).unsqueeze(1)


def periodic_heat_operator(
    initial_field: Tensor,
    diffusion: float = 0.05,
    duration: float = 0.2,
) -> Tensor:
    """Apply the exact periodic heat semigroup to [batch, height, width] fields."""
    if initial_field.ndim != 3:
        raise ValueError("initial_field must have shape [batch, height, width].")
    if min(initial_field.shape[-2:]) < 4:
        raise ValueError("Each spatial grid dimension must be at least 4.")
    if diffusion <= 0 or not math.isfinite(diffusion):
        raise ValueError("diffusion must be finite and positive.")
    if duration <= 0 or not math.isfinite(duration):
        raise ValueError("duration must be finite and positive.")
    height, width = initial_field.shape[-2:]
    frequencies_x = torch.fft.fftfreq(
        height, d=1.0 / height, device=initial_field.device
    ) * (2.0 * math.pi)
    frequencies_y = torch.fft.fftfreq(
        width, d=1.0 / width, device=initial_field.device
    ) * (2.0 * math.pi)
    wave_x, wave_y = torch.meshgrid(
        frequencies_x, frequencies_y, indexing="ij"
    )
    multiplier = torch.exp(
        -diffusion * duration * (wave_x.square() + wave_y.square())
    )
    transformed = torch.fft.fftn(initial_field, dim=(-2, -1))
    output = torch.fft.ifftn(
        transformed * multiplier, dim=(-2, -1)
    ).real
    return output
