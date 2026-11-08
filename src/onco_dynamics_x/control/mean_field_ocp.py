"""Optional CasADi optimal-control experiment for a dimensionless mean-field toy model."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ControlExperiment:
    """Output of a simulation-only optimization; control values are not doses."""

    time: np.ndarray
    cancer_proxy: np.ndarray
    immune_proxy: np.ndarray
    abstract_input: np.ndarray
    objective_value: float


def optimize_hypothetical_input(
    horizon: float = 1.0,
    intervals: int = 40,
    maximum_abstract_input: float = 0.5,
    cancer_weight: float = 1.0,
    input_penalty: float = 0.1,
) -> ControlExperiment:
    """Minimize a toy objective over an abstract control in a 2-state ODE."""
    if horizon <= 0 or intervals < 2 or maximum_abstract_input < 0:
        raise ValueError("Horizon, intervals and input bound are invalid.")
    if cancer_weight < 0 or input_penalty < 0:
        raise ValueError("Objective weights must be non-negative.")
    try:
        import casadi as ca
    except ImportError as error:
        raise ImportError(
            "Install the control extra with pip install -e .[control]."
        ) from error

    dt = horizon / intervals
    optimizer = ca.Opti()
    cancer = optimizer.variable(intervals + 1)
    immune = optimizer.variable(intervals + 1)
    control = optimizer.variable(intervals)
    optimizer.subject_to(cancer[0] == 0.35)
    optimizer.subject_to(immune[0] == 0.15)
    optimizer.subject_to(optimizer.bounded(0.0, cancer, 1.0))
    optimizer.subject_to(optimizer.bounded(0.0, immune, 2.0))
    optimizer.subject_to(
        optimizer.bounded(0.0, control, maximum_abstract_input)
    )

    objective = 0
    growth = 0.45
    immune_killing = 0.25
    abstract_control_effect = 0.35
    immune_source = 0.08
    immune_decay = 0.15

    for index in range(intervals):
        c = cancer[index]
        i = immune[index]
        u = control[index]
        cancer_rate = (
            growth * c * (1.0 - c)
            - immune_killing * i * c
            - abstract_control_effect * u * c
        )
        immune_rate = immune_source + 0.10 * c / (0.2 + c) - immune_decay * i
        optimizer.subject_to(cancer[index + 1] == c + dt * cancer_rate)
        optimizer.subject_to(immune[index + 1] == i + dt * immune_rate)
        objective += dt * (cancer_weight * c + input_penalty * u**2)
    objective += cancer_weight * cancer[-1]
    optimizer.minimize(objective)
    optimizer.solver("ipopt", {"print_time": False}, {"print_level": 0})
    solution = optimizer.solve()

    return ControlExperiment(
        time=np.linspace(0.0, horizon, intervals + 1),
        cancer_proxy=np.asarray(solution.value(cancer)).reshape(-1),
        immune_proxy=np.asarray(solution.value(immune)).reshape(-1),
        abstract_input=np.asarray(solution.value(control)).reshape(-1),
        objective_value=float(solution.value(objective)),
    )
