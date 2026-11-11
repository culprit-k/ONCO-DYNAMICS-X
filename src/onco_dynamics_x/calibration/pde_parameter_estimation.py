"""Low-dimensional PDE parameter estimation against registered field observations.

All parameters estimated here are dimensionless ModelParameters coefficients.
The observation arrays must already be registered to the model grid and time
scale. This routine does not infer image registration, physical voxel scales,
or biological equivalence of a proxy field.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from typing import Any

import numpy as np
from numpy.typing import NDArray

from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.solver import step
from onco_dynamics_x.state import TumourState

FloatArray = NDArray[np.float64]
_ALLOWED_PARAMETER_NAMES = {
    "cancer_growth",
    "baseline_cancer_loss",
    "immune_killing",
    "oxygen_consumption",
    "oxygen_supply",
    "nutrient_consumption",
    "nutrient_supply",
    "fibroblast_activation",
    "fibroblast_decay",
    "evasion_production",
    "evasion_decay",
    "vascular_adaptation",
    "vascular_return",
    "drug_decay",
}
_ALLOWED_OBSERVATION_FIELDS = set(TumourState.FIELD_NAMES)


@dataclass(frozen=True)
class PDEParameterFit:
    """Fit summary with local identifiability and covariance diagnostics."""

    parameter_estimates: dict[str, float]
    parameter_standard_errors: dict[str, float] | None
    parameter_covariance: list[list[float]] | None
    observation_fields: tuple[str, ...]
    observation_count: int
    residual_rmse_standardised: float
    jacobian_singular_values: list[float]
    jacobian_rank: int
    jacobian_condition_number: float
    locally_identifiable: bool
    optimiser_success: bool
    optimiser_message: str
    parameter_space: str = "dimensionless ModelParameters coefficients"
    uncertainty_method: str = (
        "local linearised covariance from weighted least squares; conditional "
        "on fixed initial state, observation noise scales and model structure"
    )
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _state_from_observations(
    observations: dict[str, FloatArray],
    index: int,
    template: TumourState,
) -> TumourState:
    fields = {
        name: (
            np.asarray(observations[name][index], dtype=float).copy()
            if name in observations
            else getattr(template, name).copy()
        )
        for name in TumourState.FIELD_NAMES
    }
    return TumourState(**fields)


def _simulate_observation_fields(
    initial: TumourState,
    parameters: ModelParameters,
    observation_times: FloatArray,
    dt: float,
    dx: float,
    fields: tuple[str, ...],
) -> dict[str, FloatArray]:
    ratios = np.diff(observation_times) / dt
    integer_steps = np.rint(ratios).astype(int)
    if np.any(integer_steps < 1) or not np.allclose(ratios, integer_steps, atol=1e-8):
        raise ValueError("Observation times must be positive integer multiples of dt.")
    current = initial.copy()
    output = {name: [getattr(current, name).copy()] for name in fields}
    for interval_steps in integer_steps:
        for _ in range(int(interval_steps)):
            current = step(current, parameters, dx, dt)
        for name in fields:
            output[name].append(getattr(current, name).copy())
    return {name: np.stack(values, axis=0) for name, values in output.items()}


def fit_pde_parameter_subset(
    observation_times: NDArray[np.float64] | list[float],
    observations: dict[str, NDArray[np.float64]],
    initial_state: TumourState,
    parameters: ModelParameters,
    *,
    parameter_names: tuple[str, ...] | list[str],
    parameter_bounds: dict[str, tuple[float, float]],
    observation_sd: dict[str, float],
    dx: float,
    dt: float,
    max_nfev: int = 300,
    identifiability_condition_limit: float = 1e8,
) -> PDEParameterFit:
    """Fit at most three positive reaction parameters to registered field data.

    The first observation time must be zero and all observation arrays must
    include the initial state. The initial state is held fixed. Parameters not
    named in parameter_names remain fixed at the supplied baseline. Weighted
    residuals use analyst-supplied observation standard deviations.

    This is suitable only after spatial registration, field definitions, units,
    time scaling and observation-error scales have been independently justified.
    """
    try:
        from scipy.optimize import least_squares
    except ImportError as error:
        raise ImportError(
            "Install the science extra: pip install -e '.[science]'."
        ) from error

    times = np.asarray(observation_times, dtype=float)
    names = tuple(parameter_names)
    if times.ndim != 1 or len(times) < 2 or not np.isfinite(times).all():
        raise ValueError("observation_times must be a finite 1D array with at least two points.")
    if not np.isclose(times[0], 0.0) or np.any(np.diff(times) <= 0):
        raise ValueError("Observation times must start at zero and be strictly increasing.")
    if dx <= 0 or dt <= 0 or max_nfev < 1:
        raise ValueError("dx, dt and max_nfev must be positive.")
    if not 1 <= len(names) <= 3 or len(set(names)) != len(names):
        raise ValueError("Fit one to three unique parameters per calibration stage.")
    unknown = set(names) - _ALLOWED_PARAMETER_NAMES
    if unknown:
        raise ValueError(f"Unsupported calibration parameter(s): {', '.join(sorted(unknown))}.")
    if set(names) != set(parameter_bounds):
        raise ValueError("parameter_bounds must define exactly the fitted parameters.")
    if not observations:
        raise ValueError("At least one observed field is required.")
    fields = tuple(observations)
    if set(fields) - _ALLOWED_OBSERVATION_FIELDS:
        raise ValueError("Observations contain an unknown TumourState field.")
    if set(fields) != set(observation_sd):
        raise ValueError("Provide one observation_sd for every observed field.")

    arrays: dict[str, FloatArray] = {}
    for field in fields:
        array = np.asarray(observations[field], dtype=float)
        if array.shape != (len(times), *initial_state.shape):
            raise ValueError(
                f"{field} observations must have shape (time, *grid) matching initial_state."
            )
        if not np.isfinite(array).all():
            raise ValueError(f"{field} observations must be finite.")
        scale = float(observation_sd[field])
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError(f"observation_sd for {field} must be finite and positive.")
        if not np.allclose(array[0], getattr(initial_state, field), rtol=1e-6, atol=1e-8):
            raise ValueError(
                f"{field} observation at time zero must match initial_state."
            )
        arrays[field] = array

    lower, upper, start = [], [], []
    for name in names:
        if name not in parameter_bounds:
            raise ValueError(f"Missing bounds for {name}.")
        low, high = map(float, parameter_bounds[name])
        baseline = float(getattr(parameters, name))
        if not np.isfinite([low, high, baseline]).all() or low <= 0 or high <= low:
            raise ValueError(f"Bounds for {name} must be finite, positive and ordered.")
        if not low <= baseline <= high:
            raise ValueError(f"Baseline value for {name} must lie within its bounds.")
        lower.append(np.log(low))
        upper.append(np.log(high))
        start.append(np.log(baseline))

    observed_vector = np.concatenate(
        [arrays[field].reshape(-1) / observation_sd[field] for field in fields]
    )

    def predict(log_parameters: FloatArray) -> FloatArray:
        changes = {
            name: float(np.exp(value))
            for name, value in zip(names, log_parameters, strict=True)
        }
        candidate = replace(parameters, **changes)
        predicted = _simulate_observation_fields(
            initial_state, candidate, times, dt, dx, fields
        )
        return np.concatenate(
            [
                predicted[field].reshape(-1) / observation_sd[field]
                for field in fields
            ]
        )

    fit = least_squares(
        lambda theta: predict(theta) - observed_vector,
        x0=np.asarray(start),
        bounds=(np.asarray(lower), np.asarray(upper)),
        method="trf",
        max_nfev=max_nfev,
        x_scale="jac",
    )
    residuals = fit.fun
    singular_values = np.linalg.svd(fit.jac, compute_uv=False)
    tolerance = (
        np.finfo(float).eps * max(fit.jac.shape) * singular_values[0]
        if singular_values.size
        else 0.0
    )
    rank = int(np.sum(singular_values > tolerance))
    condition = (
        float(singular_values[0] / singular_values[-1])
        if singular_values.size and singular_values[-1] > tolerance
        else float("inf")
    )
    dof = max(len(residuals) - len(names), 1)
    residual_variance = float(np.dot(residuals, residuals) / dof)
    covariance_log = None
    standard_errors = None
    warnings: list[str] = []
    identifiable = rank == len(names) and condition < identifiability_condition_limit
    if identifiable:
        covariance_log = residual_variance * np.linalg.inv(fit.jac.T @ fit.jac)
        estimates = {
            name: float(np.exp(value))
            for name, value in zip(names, fit.x, strict=True)
        }
        # Delta-method conversion from log-parameter to natural parameter scale.
        scale_matrix = np.diag([estimates[name] for name in names])
        covariance = scale_matrix @ covariance_log @ scale_matrix
        standard_errors = {
            name: float(np.sqrt(max(covariance[index, index], 0.0)))
            for index, name in enumerate(names)
        }
        covariance_result = covariance.tolist()
    else:
        estimates = {
            name: float(np.exp(value))
            for name, value in zip(names, fit.x, strict=True)
        }
        covariance_result = None
        warnings.append(
            "The fitted parameter subset is locally non-identifiable or ill-conditioned; "
            "do not interpret point estimates as separately determined biological values."
        )
    if not fit.success:
        warnings.append("The nonlinear optimiser did not report convergence.")
    if np.any(np.isclose(fit.x, lower, atol=1e-5)) or np.any(
        np.isclose(fit.x, upper, atol=1e-5)
    ):
        warnings.append("At least one estimate is at or near a parameter bound.")
    return PDEParameterFit(
        parameter_estimates=estimates,
        parameter_standard_errors=standard_errors,
        parameter_covariance=covariance_result,
        observation_fields=fields,
        observation_count=int(observed_vector.size),
        residual_rmse_standardised=float(np.sqrt(np.mean(residuals**2))),
        jacobian_singular_values=[float(value) for value in singular_values],
        jacobian_rank=rank,
        jacobian_condition_number=condition,
        locally_identifiable=identifiable,
        optimiser_success=bool(fit.success),
        optimiser_message=str(fit.message),
        warnings=tuple(warnings),
    )
