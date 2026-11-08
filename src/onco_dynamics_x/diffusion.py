"""Finite-difference spatial operators with no-flux edge conditions."""
from __future__ import annotations
import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]

def laplacian_no_flux(field: FloatArray, dx: float) -> FloatArray:
    """Return the 2D five-point Laplacian with edge-value (zero-flux) padding."""
    if field.ndim != 2:
        raise ValueError("field must be a 2D array.")
    if dx <= 0:
        raise ValueError("dx must be positive.")
    padded = np.pad(field, pad_width=1, mode="edge")
    centre = padded[1:-1, 1:-1]
    north, south = padded[:-2, 1:-1], padded[2:, 1:-1]
    west, east = padded[1:-1, :-2], padded[1:-1, 2:]
    return (north + south + west + east - 4.0 * centre) / dx**2

def diffusion_term(field: FloatArray, coefficient: float, dx: float) -> FloatArray:
    if coefficient < 0:
        raise ValueError("Diffusion coefficient must be non-negative.")
    return coefficient * laplacian_no_flux(field, dx)
