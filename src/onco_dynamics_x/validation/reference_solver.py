"""Independent adaptive method-of-lines reference for the reaction-diffusion system.

This solver shares the declared spatial finite-difference operator with the
explicit solver but uses SciPy's adaptive ODE integration, providing an
independent time-integration check on small grids. It is a numerical reference,
not a biological ground truth.
"""
from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

from onco_dynamics_x.diffusion import diffusion_term
from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.reactions import reaction_terms
from onco_dynamics_x.state import TumourState
from onco_dynamics_x.validation.diagnostics import relative_l2_error


def simulate_adaptive_reference(
    initial: TumourState,
    parameters: ModelParameters,
    dx: float,
    final_time: float,
    *,
    method: str = "DOP853",
    rtol: float = 1e-8,
    atol: float = 1e-10,
    max_step: float = np.inf,
    infusion_schedule: Callable[[float], float] | None = None,
) -> tuple[TumourState, dict[str, Any]]:
    """Integrate the semi-discrete 2D/3D model with an adaptive SciPy solver."""
    try:
        from scipy.integrate import solve_ivp
    except ImportError as exc:
        raise RuntimeError(
            "Install the science dependencies with pip install -e '.[science]'."
        ) from exc

    initial.validate()
    parameters.validate()
    if dx <= 0 or not np.isfinite(dx):
        raise ValueError("dx must be finite and positive.")
    if final_time <= 0 or not np.isfinite(final_time):
        raise ValueError("final_time must be finite and positive.")
    if rtol <= 0 or atol <= 0:
        raise ValueError("rtol and atol must be positive.")
    if max_step <= 0:
        raise ValueError("max_step must be positive.")

    field_names = initial.FIELD_NAMES
    field_shape = initial.shape
    field_size = initial.cancer.size
    initial_vector = np.concatenate(
        [initial.as_dict()[name].ravel() for name in field_names]
    )
    diffusion = parameters.diffusion_coefficients

    def unpack(vector: np.ndarray) -> TumourState:
        fields = {}
        for index, name in enumerate(field_names):
            start = index * field_size
            stop = start + field_size
            fields[name] = vector[start:stop].reshape(field_shape)
        return TumourState(**fields)

    def right_hand_side(time: float, vector: np.ndarray) -> np.ndarray:
        state = unpack(vector)
        input_rate = (
            0.0 if infusion_schedule is None else float(infusion_schedule(time))
        )
        if not np.isfinite(input_rate) or input_rate < 0:
            raise ValueError("infusion_schedule must return finite non-negative values.")
        reactions = reaction_terms(state, parameters, infusion_rate=input_rate)
        derivatives = []
        for name, field in state.as_dict().items():
            derivatives.append(
                (diffusion_term(field, diffusion[name], dx) + reactions[name]).ravel()
            )
        return np.concatenate(derivatives)

    solution = solve_ivp(
        right_hand_side,
        (0.0, float(final_time)),
        initial_vector,
        method=method,
        rtol=rtol,
        atol=atol,
        max_step=max_step,
        t_eval=[float(final_time)],
    )
    if not solution.success:
        raise RuntimeError(f"Adaptive reference integration failed: {solution.message}")
    final_state = unpack(solution.y[:, -1])
    final_state.validate()
    diagnostics = {
        "method": method,
        "success": bool(solution.success),
        "message": str(solution.message),
        "final_time": float(final_time),
        "rtol": float(rtol),
        "atol": float(atol),
        "max_step": float(max_step),
        "function_evaluations": int(solution.nfev),
        "jacobian_evaluations": int(solution.njev),
        "lu_decompositions": int(solution.nlu),
        "spatial_dimensions": int(initial.cancer.ndim),
        "grid_shape": list(initial.shape),
        "post_integration_clipping_applied": False,
    }
    return final_state, diagnostics


def compare_explicit_euler_to_reference(
    initial: TumourState,
    parameters: ModelParameters,
    dx: float,
    dt: float,
    steps: int,
    *,
    method: str = "DOP853",
    rtol: float = 1e-8,
    atol: float = 1e-10,
    infusion_schedule: Callable[[float], float] | None = None,
) -> dict[str, Any]:
    """Compare explicit Euler and adaptive ODE integration at equal final time."""
    if steps < 1 or dt <= 0:
        raise ValueError("steps and dt must be positive.")

    from onco_dynamics_x.solver import simulate

    _, snapshots = simulate(
        initial,
        parameters,
        dx=dx,
        dt=dt,
        steps=steps,
        infusion_schedule=infusion_schedule,
        save_every=steps,
    )
    reference, diagnostics = simulate_adaptive_reference(
        initial,
        parameters,
        dx=dx,
        final_time=dt * steps,
        method=method,
        rtol=rtol,
        atol=atol,
        infusion_schedule=infusion_schedule,
    )
    explicit = snapshots[-1].as_dict()
    adaptive = reference.as_dict()
    errors = {
        name: relative_l2_error(adaptive[name], explicit[name])
        for name in initial.FIELD_NAMES
    }
    return {
        "comparison": "explicit_euler_vs_adaptive_method_of_lines",
        "final_time": float(dt * steps),
        "explicit_dt": float(dt),
        "explicit_steps": int(steps),
        "relative_l2_error_by_field": errors,
        "mean_relative_l2_error": float(np.mean(list(errors.values()))),
        "maximum_relative_l2_error": float(max(errors.values())),
        "adaptive_reference": diagnostics,
        "interpretation": (
            "Both solvers use the same spatial finite-difference operator and "
            "reaction terms; this comparison checks time integration, not spatial "
            "discretization or biological correctness."
        ),
    }
