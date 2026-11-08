"""Explicit Euler reaction-diffusion solver for the 2D prototype."""
from __future__ import annotations
from collections.abc import Callable
import numpy as np
from onco_dynamics_x.diffusion import diffusion_term
from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.reactions import reaction_terms
from onco_dynamics_x.state import TumourState

def diffusion_timestep_limit(parameters: ModelParameters, dx: float) -> float:
    """2D explicit diffusion limit dx^2 / (4 * max(D))."""
    if dx <= 0:
        raise ValueError("dx must be positive.")
    maximum = max(parameters.diffusion_coefficients.values())
    return float("inf") if maximum == 0 else dx**2 / (4.0 * maximum)

def _bound_state(state: TumourState, parameters: ModelParameters) -> TumourState:
    p = parameters
    state.cancer = np.clip(state.cancer, 0.0, p.cancer_capacity)
    state.immune = np.clip(state.immune, 0.0, p.maximum_field)
    state.fibroblast = np.clip(state.fibroblast, 0.0, p.maximum_field)
    state.oxygen = np.clip(state.oxygen, 0.0, p.maximum_field)
    state.nutrient = np.clip(state.nutrient, 0.0, p.maximum_field)
    state.drug = np.clip(state.drug, 0.0, p.maximum_drug)
    state.evasion = np.clip(state.evasion, 0.0, p.maximum_field)
    state.vasculature = np.clip(state.vasculature, 0.1, p.maximum_field)
    return state

def step(
    state: TumourState, parameters: ModelParameters, dx: float, dt: float,
    infusion_rate: float = 0.0, *, enforce_diffusion_stability: bool = True
) -> TumourState:
    """Advance one explicit Euler step.

    The stability guard covers diffusion only; reaction stiffness requires
    additional timestep-refinement and convergence analysis.
    """
    if dx <= 0 or dt <= 0:
        raise ValueError("dx and dt must be positive.")
    parameters.validate()
    state.validate()
    limit = diffusion_timestep_limit(parameters, dx)
    if enforce_diffusion_stability and dt > limit * (1.0 + 1e-12):
        raise ValueError(
            f"dt={dt:g} exceeds the explicit diffusion limit {limit:g}; reduce dt or refine the grid."
        )
    reactions = reaction_terms(state, parameters, infusion_rate)
    next_fields = {}
    for name, field in state.as_dict().items():
        diffusion = diffusion_term(field, parameters.diffusion_coefficients[name], dx)
        next_fields[name] = field + dt * (diffusion + reactions[name])
    return _bound_state(TumourState(**next_fields), parameters)

def simulate(
    initial: TumourState, parameters: ModelParameters, dx: float, dt: float, steps: int,
    infusion_schedule: Callable[[float], float] | None = None, save_every: int = 10
) -> tuple[np.ndarray, list[TumourState]]:
    """Run a deterministic simulation and retain regularly spaced state snapshots."""
    if steps < 1:
        raise ValueError("steps must be at least 1.")
    if save_every < 1:
        raise ValueError("save_every must be at least 1.")
    current = initial.copy()
    current.validate()
    snapshots, times = [current.copy()], [0.0]
    for index in range(1, steps + 1):
        time = (index - 1) * dt
        input_rate = 0.0 if infusion_schedule is None else float(infusion_schedule(time))
        current = step(current, parameters, dx, dt, input_rate)
        if index % save_every == 0 or index == steps:
            snapshots.append(current.copy())
            times.append(index * dt)
    return np.asarray(times, dtype=float), snapshots
