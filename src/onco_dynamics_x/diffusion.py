"""Finite-difference spatial operators with no-flux edge conditions."""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


def laplacian_no_flux(field: FloatArray, dx: float) -> FloatArray:
    """Return the 2D or 3D nearest-neighbour Laplacian with edge-value padding."""
    if field.ndim not in (2, 3):
        raise ValueError("field must be a 2D or 3D array.")
    if dx <= 0:
        raise ValueError("dx must be positive.")
    padding = [(1, 1)] * field.ndim
    padded = np.pad(field, pad_width=padding, mode="edge")
    centre_slice = (slice(1, -1),) * field.ndim
    centre = padded[centre_slice]
    result = -2.0 * field.ndim * centre
    for axis in range(field.ndim):
        lower = [slice(1, -1)] * field.ndim
        upper = [slice(1, -1)] * field.ndim
        lower[axis] = slice(0, -2)
        upper[axis] = slice(2, None)
        result = result + padded[tuple(lower)] + padded[tuple(upper)]
    return result / dx**2


def diffusion_term(field: FloatArray, coefficient: float, dx: float) -> FloatArray:
    """Compute D times the spatial Laplacian."""
    if coefficient < 0:
        raise ValueError("Diffusion coefficient must be non-negative.")
    return coefficient * laplacian_no_flux(field, dx)
