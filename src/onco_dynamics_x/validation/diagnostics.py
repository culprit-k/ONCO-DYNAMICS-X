"""Numerical refinement diagnostics for the reference PDE solver."""
from __future__ import annotations

import numpy as np

from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.solver import simulate
from onco_dynamics_x.state import TumourState


def relative_l2_error(reference: np.ndarray, candidate: np.ndarray) -> float:
    """Return relative L2 error with a safe denominator for zero reference arrays."""
    if reference.shape != candidate.shape:
        raise ValueError("Reference and candidate arrays must have the same shape.")
    if not np.isfinite(reference).all() or not np.isfinite(candidate).all():
        raise ValueError("Arrays must contain only finite values.")
    denominator = float(np.linalg.norm(reference.ravel()))
    numerator = float(np.linalg.norm((candidate - reference).ravel()))
    return numerator if denominator == 0.0 else numerator / denominator


def timestep_refinement(
    initial: TumourState,
    parameters: ModelParameters,
    dx: float,
    dt: float,
    steps: int,
) -> dict[str, float]:
    """Compare dt to dt/2 over the same final time and return per-field errors."""
    if steps < 1 or dt <= 0:
        raise ValueError("steps and dt must be positive.")
    _, coarse_snapshots = simulate(initial, parameters, dx, dt, steps, save_every=steps)
    _, fine_snapshots = simulate(initial, parameters, dx, dt / 2.0, steps * 2, save_every=steps * 2)
    coarse = coarse_snapshots[-1].as_dict()
    fine = fine_snapshots[-1].as_dict()
    return {
        name: relative_l2_error(coarse[name], fine[name])
        for name in initial.FIELD_NAMES
    }
