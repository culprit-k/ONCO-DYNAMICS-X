"""Numerical timestep and grid-refinement diagnostics for the reference PDE solver."""
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
    _, fine_snapshots = simulate(
        initial, parameters, dx, dt / 2.0, steps * 2, save_every=steps * 2
    )
    coarse = coarse_snapshots[-1].as_dict()
    fine = fine_snapshots[-1].as_dict()
    return {
        name: relative_l2_error(coarse[name], fine[name])
        for name in initial.FIELD_NAMES
    }


def grid_refinement(
    coarse_initial: TumourState,
    fine_initial: TumourState,
    parameters: ModelParameters,
    coarse_dx: float,
    coarse_dt: float,
    coarse_steps: int,
    refinement_ratio: int = 2,
) -> dict[str, float]:
    """Compare nested coarse/fine grids at equal final time.

    The fine grid must include the coarse nodes exactly: n_fine=(n_coarse-1)*r+1.
    The fine timestep is reduced by r^2 and the step count increased by r^2.
    Fine outputs are sampled at coarse-grid locations before relative-L2 comparison.
    """
    if refinement_ratio < 2:
        raise ValueError("refinement_ratio must be at least 2.")
    if coarse_steps < 1 or coarse_dt <= 0 or coarse_dx <= 0:
        raise ValueError("coarse_steps, coarse_dt and coarse_dx must be positive.")
    coarse_initial.validate()
    fine_initial.validate()
    if coarse_initial.cancer.ndim != fine_initial.cancer.ndim:
        raise ValueError("Coarse and fine grids must have the same spatial dimension.")
    dimensions = coarse_initial.cancer.ndim
    expected_size = (coarse_initial.shape[0] - 1) * refinement_ratio + 1
    expected_shape = (expected_size,) * dimensions
    if fine_initial.shape != expected_shape:
        raise ValueError(
            f"Fine grid shape must be {expected_shape} for refinement ratio "
            f"{refinement_ratio}."
        )
    fine_dx = coarse_dx / refinement_ratio
    fine_dt = coarse_dt / (refinement_ratio**2)
    fine_steps = coarse_steps * refinement_ratio**2

    _, coarse_snapshots = simulate(
        coarse_initial, parameters, coarse_dx, coarse_dt, coarse_steps,
        save_every=coarse_steps,
    )
    _, fine_snapshots = simulate(
        fine_initial, parameters, fine_dx, fine_dt, fine_steps,
        save_every=fine_steps,
    )
    coarse = coarse_snapshots[-1].as_dict()
    fine = fine_snapshots[-1].as_dict()
    slices = (slice(None, None, refinement_ratio),) * dimensions
    return {
        name: relative_l2_error(coarse[name], fine[name][slices])
        for name in coarse_initial.FIELD_NAMES
    }
