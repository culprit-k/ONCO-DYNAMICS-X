"""Local finite-difference sensitivity analysis for scalar simulation summaries."""
from __future__ import annotations

from dataclasses import replace
from typing import Callable

import numpy as np

from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.solver import simulate
from onco_dynamics_x.state import TumourState


def final_mean_cancer(
    parameters: ModelParameters,
    initial: TumourState,
    dx: float,
    dt: float,
    steps: int,
) -> float:
    """Evaluate mean cancer proxy at final time for one simulation."""
    _, snapshots = simulate(initial, parameters, dx, dt, steps, save_every=steps)
    return float(np.mean(snapshots[-1].cancer))


def local_parameter_sensitivity(
    parameter_name: str,
    parameters: ModelParameters,
    initial: TumourState,
    dx: float,
    dt: float,
    steps: int,
    relative_step: float = 1e-3,
    metric: Callable[[ModelParameters, TumourState, float, float, int], float] = final_mean_cancer,
) -> dict[str, float]:
    """Estimate a local normalized sensitivity by a central finite difference.

    The returned elasticity is (parameter / metric) * d(metric)/d(parameter).
    It is a numerical property of the chosen model/configuration, not a biological
    causal effect or a population-level estimate.
    """
    if parameter_name not in parameters.__dataclass_fields__:
        raise ValueError(f"Unknown model parameter: {parameter_name}")
    if relative_step <= 0:
        raise ValueError("relative_step must be positive.")
    value = float(getattr(parameters, parameter_name))
    if value <= 0:
        raise ValueError("Sensitivity requires a positive baseline parameter.")
    increment = value * relative_step
    low = max(value - increment, value * 0.5)
    high = value + increment
    low_parameters = replace(parameters, **{parameter_name: low})
    high_parameters = replace(parameters, **{parameter_name: high})
    base_metric = metric(parameters, initial, dx, dt, steps)
    low_metric = metric(low_parameters, initial, dx, dt, steps)
    high_metric = metric(high_parameters, initial, dx, dt, steps)
    derivative = (high_metric - low_metric) / (high - low)
    elasticity = 0.0 if base_metric == 0 else value * derivative / base_metric
    return {
        "parameter_value": value,
        "metric_baseline": base_metric,
        "derivative": derivative,
        "elasticity": elasticity,
        "finite_difference_step": high - low,
    }
