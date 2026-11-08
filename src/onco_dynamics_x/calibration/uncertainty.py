"""Monte Carlo parameter-range uncertainty for deterministic simulations."""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from numpy.typing import NDArray

from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.solver import simulate
from onco_dynamics_x.state import TumourState

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class ParameterUncertaintyResult:
    """Distribution of a simulation summary under stated parameter ranges."""

    parameter_samples: dict[str, FloatArray]
    final_cancer_mean: FloatArray
    mean: float
    median: float
    lower_95: float
    upper_95: float


def monte_carlo_parameter_uncertainty(
    initial: TumourState,
    parameters: ModelParameters,
    dx: float,
    dt: float,
    steps: int,
    parameter_ranges: dict[str, tuple[float, float]],
    samples: int = 100,
    seed: int = 42,
) -> ParameterUncertaintyResult:
    """Propagate user-specified uniform parameter ranges through the PDE model.

    This is range-based Monte Carlo uncertainty, not a Bayesian posterior. Bounds
    must be justified externally and must not be interpreted as clinical ranges.
    """
    if samples < 2 or steps < 1 or dx <= 0 or dt <= 0:
        raise ValueError("samples, steps, dx and dt must be positive (samples >= 2).")
    if not parameter_ranges:
        raise ValueError("At least one parameter range is required.")
    allowed = set(parameters.__dataclass_fields__)
    unknown = set(parameter_ranges) - allowed
    if unknown:
        raise ValueError(f"Unknown parameter(s): {', '.join(sorted(unknown))}.")

    rng = np.random.default_rng(seed)
    parameter_samples: dict[str, FloatArray] = {}
    for name, bounds in parameter_ranges.items():
        if len(bounds) != 2:
            raise ValueError(f"Range for {name} must contain (lower, upper).")
        lower, upper = map(float, bounds)
        if not np.isfinite([lower, upper]).all() or lower > upper:
            raise ValueError(f"Range for {name} must be finite with lower <= upper.")
        if lower < 0:
            raise ValueError(f"Range for {name} cannot include negative values.")
        parameter_samples[name] = rng.uniform(lower, upper, samples)

    outcomes = np.empty(samples, dtype=float)
    for index in range(samples):
        changes = {
            name: float(values[index]) for name, values in parameter_samples.items()
        }
        current_parameters = replace(parameters, **changes)
        current_parameters.validate()
        _, snapshots = simulate(
            initial, current_parameters, dx, dt, steps, save_every=steps
        )
        outcomes[index] = float(np.mean(snapshots[-1].cancer))

    return ParameterUncertaintyResult(
        parameter_samples=parameter_samples,
        final_cancer_mean=outcomes,
        mean=float(np.mean(outcomes)),
        median=float(np.median(outcomes)),
        lower_95=float(np.quantile(outcomes, 0.025)),
        upper_95=float(np.quantile(outcomes, 0.975)),
    )
