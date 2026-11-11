"""Analytic synthetic operators used to benchmark optional AI architectures.

These operators provide known mathematical ground truth for software tests.
They are not breast-cancer models or patient-derived training data.
"""
from __future__ import annotations

import math

import numpy as np
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



def synthetic_spatial_neighbour_graph(
    num_nodes: int = 64,
    radius: float = 0.4,
    seed: int = 42,
) -> tuple[Tensor, Tensor, Tensor]:
    """Create a synthetic cell graph with known local immune-neighbour fractions.

    Node features contain x/y coordinates and one-hot cancer/immune/fibroblast
    phenotypes. The scalar target is the fraction of immune-type neighbours
    within the fixed Euclidean radius. This is a software benchmark, not a
    model of biological cell recruitment or a claim about a tumour microenvironment.
    """
    if num_nodes < 4:
        raise ValueError("num_nodes must be at least 4.")
    if not math.isfinite(radius) or radius <= 0:
        raise ValueError("radius must be finite and positive.")

    rng = np.random.default_rng(seed)
    coordinates = rng.uniform(-1.0, 1.0, size=(num_nodes, 2))
    phenotypes = rng.integers(0, 3, size=num_nodes)
    one_hot = np.eye(3, dtype=float)[phenotypes]
    features = np.column_stack([coordinates, one_hot])

    differences = coordinates[:, None, :] - coordinates[None, :, :]
    distance_squared = np.sum(differences**2, axis=2)
    adjacency = (distance_squared <= radius**2) & ~np.eye(
        num_nodes, dtype=bool
    )
    source, target = np.nonzero(adjacency)
    immune_indicator = (phenotypes == 1).astype(float)
    immune_neighbour_count = np.bincount(
        target,
        weights=immune_indicator[source],
        minlength=num_nodes,
    )
    neighbour_count = np.bincount(target, minlength=num_nodes)
    neighbour_fraction = np.divide(
        immune_neighbour_count,
        neighbour_count,
        out=np.zeros(num_nodes, dtype=float),
        where=neighbour_count > 0,
    )
    edge_index = np.vstack([source, target]).astype(np.int64, copy=False)
    return (
        torch.as_tensor(features, dtype=torch.float32),
        torch.as_tensor(edge_index, dtype=torch.long),
        torch.as_tensor(neighbour_fraction[:, None], dtype=torch.float32),
    )
