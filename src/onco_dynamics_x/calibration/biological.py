"""Estimation of effective tumour-volume kinetics from longitudinal observations.

This module fits scalar untreated-growth curves to repeated physical tumour
volume measurements. Its rates are effective phenomenological quantities in
1/day, not direct estimates of the dimensionless PDE coefficients in
ModelParameters. Therapy-exposed intervals must not be fitted as untreated
growth data.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Literal

import numpy as np
from numpy.typing import NDArray

GrowthModel = Literal["exponential", "logistic", "gompertz"]
FloatArray = NDArray[np.float64]
_ALLOWED_MODELS: tuple[GrowthModel, ...] = (
    "exponential",
    "logistic",
    "gompertz",
)


@dataclass(frozen=True)
class GrowthCurveFit:
    """Single-series fit and local identifiability diagnostics."""

    model: str
    parameters: dict[str, float]
    time_unit: str
    volume_unit: str
    n_observations: int
    residual_rmse_log: float
    residual_sd_log: float
    aicc: float | None
    bic: float
    jacobian_singular_values: list[float]
    jacobian_condition_number: float
    locally_identifiable: bool
    carrying_capacity_weakly_identified: bool
    bootstrap_95ci: dict[str, list[float]]
    warnings: list[str]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _predict(
    model: GrowthModel,
    time_days: FloatArray,
    initial_volume: float,
    theta: FloatArray,
    maximum_observed: float,
) -> tuple[FloatArray, dict[str, float]]:
    rate = float(np.exp(theta[0]))
    if model == "exponential":
        predicted = initial_volume * np.exp(rate * time_days)
        return predicted, {"growth_rate_per_day": rate}
    capacity = maximum_observed * (1.0 + float(np.exp(theta[1])))
    if model == "logistic":
        ratio = capacity / initial_volume - 1.0
        predicted = capacity / (1.0 + ratio * np.exp(-rate * time_days))
        parameter_name = "logistic_rate_per_day"
    else:
        predicted = capacity * np.exp(
            np.log(initial_volume / capacity) * np.exp(-rate * time_days)
        )
        parameter_name = "gompertz_rate_per_day"
    return predicted, {
        parameter_name: rate,
        "carrying_capacity_mm3": capacity,
    }


def _fit_one_model(
    model: GrowthModel,
    time_days: FloatArray,
    volume_mm3: FloatArray,
) -> dict[str, Any]:
    try:
        from scipy.optimize import least_squares
    except ImportError as error:
        raise ImportError(
            "Install the science extra: pip install -e '.[science]'."
        ) from error

    initial = float(volume_mm3[0])
    maximum = float(np.max(volume_mm3))
    log_slope = float(
        np.polyfit(time_days, np.log(volume_mm3), deg=1)[0]
    )
    rate_start = float(np.clip(log_slope, 1e-4, 1.0))
    log_rate_bounds = (np.log(1e-5), np.log(5.0))

    if model == "exponential":
        starts = [np.array([np.log(rate_start)], dtype=float)]
        lower = np.array([log_rate_bounds[0]])
        upper = np.array([log_rate_bounds[1]])
    else:
        starts = [
            np.array([np.log(rate_start), np.log(margin)], dtype=float)
            for margin in (0.1, 0.5, 2.0, 10.0, 100.0)
        ]
        lower = np.array([log_rate_bounds[0], -8.0])
        upper = np.array([log_rate_bounds[1], 10.0])

    best = None
    for start in starts:
        fit = least_squares(
            lambda theta: (
                np.log(
                    _predict(model, time_days, initial, theta, maximum)[0]
                )
                - np.log(volume_mm3)
            ),
            x0=np.clip(start, lower + 1e-9, upper - 1e-9),
            bounds=(lower, upper),
            method="trf",
            max_nfev=10_000,
            ftol=1e-11,
            xtol=1e-11,
            gtol=1e-11,
        )
        if best is None or np.dot(fit.fun, fit.fun) < np.dot(best.fun, best.fun):
            best = fit
    assert best is not None

    predicted, parameters = _predict(
        model, time_days, initial, best.x, maximum
    )
    residuals = np.log(volume_mm3) - np.log(predicted)
    sse = max(float(np.dot(residuals, residuals)), 1e-15)
    n = len(volume_mm3)
    parameter_count = len(best.x)
    likelihood_parameter_count = parameter_count + 1  # include residual scale
    aic = n * np.log(sse / n) + 2.0 * likelihood_parameter_count
    aicc = (
        float(
            aic
            + 2.0 * likelihood_parameter_count
            * (likelihood_parameter_count + 1)
            / (n - likelihood_parameter_count - 1)
        )
        if n > likelihood_parameter_count + 1
        else None
    )
    bic = float(
        n * np.log(sse / n)
        + likelihood_parameter_count * np.log(n)
    )
    singular_values = np.linalg.svd(best.jac, compute_uv=False)
    if singular_values.size and singular_values[-1] > 0:
        condition = float(singular_values[0] / singular_values[-1])
    else:
        condition = float("inf")
    capacity = parameters.get("carrying_capacity_mm3")
    capacity_weak = bool(
        capacity is not None
        and (
            capacity > 100.0 * maximum
            or best.x[1] > 8.0
            or condition > 1e6
        )
    )
    warnings: list[str] = []
    if capacity_weak:
        warnings.append(
            "Carrying capacity is weakly identified; the observed window may "
            "not contain enough saturation information."
        )
    if condition > 1e6:
        warnings.append(
            "The fitted Jacobian is ill-conditioned; parameters may be confounded."
        )
    if np.any(np.isclose(best.x, lower, atol=1e-5)) or np.any(
        np.isclose(best.x, upper, atol=1e-5)
    ):
        warnings.append("At least one fitted parameter reached an optimisation bound.")

    parameters["initial_volume_mm3"] = initial
    return {
        "model": model,
        "parameters": parameters,
        "predicted": predicted,
        "residuals": residuals,
        "rmse_log": float(np.sqrt(np.mean(residuals**2))),
        "residual_sd_log": float(np.sqrt(sse / max(n - parameter_count, 1))),
        "aicc": aicc,
        "bic": bic,
        "singular_values": singular_values,
        "condition": condition,
        "locally_identifiable": bool(
            condition < 1e6 and n > parameter_count + 1
        ),
        "capacity_weak": capacity_weak,
        "warnings": warnings,
        "theta": best.x.copy(),
    }


def fit_untreated_growth_curve(
    time_days: NDArray[np.float64] | list[float],
    volume_mm3: NDArray[np.float64] | list[float],
    *,
    n_bootstrap: int = 200,
    seed: int = 42,
) -> GrowthCurveFit:
    """Compare exponential, logistic and Gompertz growth models.

    Inputs must be serial measurements from one untreated tumour, in days and
    mm^3. At least five unique time points are required. Residuals are fitted on
    log volume, corresponding to a multiplicative/log-normal measurement-error
    assumption. The selected model minimises AICc; uncertainty is estimated by
    a parametric residual bootstrap conditional on the selected model.
    """
    time = np.asarray(time_days, dtype=float)
    volume = np.asarray(volume_mm3, dtype=float)
    if time.ndim != 1 or volume.ndim != 1 or time.shape != volume.shape:
        raise ValueError("time_days and volume_mm3 must be aligned one-dimensional arrays.")
    if len(time) < 5:
        raise ValueError("At least five longitudinal observations are required.")
    if not np.isfinite(time).all() or not np.isfinite(volume).all():
        raise ValueError("Time and volume observations must be finite.")
    if np.any(volume <= 0):
        raise ValueError("Tumour volumes must be strictly positive.")
    if n_bootstrap < 0:
        raise ValueError("n_bootstrap cannot be negative.")
    order = np.argsort(time)
    time = time[order]
    volume = volume[order]
    if np.any(np.diff(time) <= 0):
        raise ValueError("Time points must be unique and strictly increasing.")
    if np.median(np.diff(np.log(volume)) / np.diff(time)) <= 0:
        raise ValueError(
            "This estimator is for untreated net growth; observed volumes do not "
            "show a positive median growth slope. Separate treatment effects first."
        )
    relative_time = time - time[0]

    fits = [_fit_one_model(model, relative_time, volume) for model in _ALLOWED_MODELS]
    finite_aicc = [fit for fit in fits if fit["aicc"] is not None]
    if not finite_aicc:
        raise ValueError("Insufficient observations to compare growth models using AICc.")
    selected = min(finite_aicc, key=lambda item: item["aicc"])

    bootstrap: dict[str, list[float]] = {}
    if n_bootstrap:
        rng = np.random.default_rng(seed)
        sampled: dict[str, list[float]] = {
            name: [] for name in selected["parameters"]
        }
        expected = selected["predicted"]
        residuals = selected["residuals"]
        centred = residuals - np.mean(residuals)
        for _ in range(n_bootstrap):
            synthetic = np.exp(
                np.log(expected) + rng.choice(centred, size=len(centred), replace=True)
            )
            try:
                replicate = _fit_one_model(
                    selected["model"], relative_time, synthetic
                )
            except (ValueError, FloatingPointError):
                continue
            for name in sampled:
                value = replicate["parameters"].get(name)
                if value is not None and np.isfinite(value):
                    sampled[name].append(float(value))
        for name, values in sampled.items():
            if len(values) >= max(30, n_bootstrap // 2):
                bootstrap[name] = [
                    float(np.quantile(values, 0.025)),
                    float(np.quantile(values, 0.975)),
                ]

    warnings = list(selected["warnings"])
    if selected["model"] != "exponential":
        warnings.append(
            "AICc preference for a saturating curve is not proof that the "
            "carrying-capacity parameter is biologically stable."
        )
    if not selected["locally_identifiable"]:
        warnings.append(
            "Local identifiability criteria were not met; do not use this fit "
            "as a fixed biological parameter."
        )
    return GrowthCurveFit(
        model=selected["model"],
        parameters=selected["parameters"],
        time_unit="day",
        volume_unit="mm^3",
        n_observations=int(len(time)),
        residual_rmse_log=selected["rmse_log"],
        residual_sd_log=selected["residual_sd_log"],
        aicc=float(selected["aicc"]),
        bic=float(selected["bic"]),
        jacobian_singular_values=[
            float(value) for value in selected["singular_values"]
        ],
        jacobian_condition_number=float(selected["condition"]),
        locally_identifiable=bool(selected["locally_identifiable"]),
        carrying_capacity_weakly_identified=bool(selected["capacity_weak"]),
        bootstrap_95ci=bootstrap,
        warnings=warnings,
    )


def physical_rate_to_dimensionless_rate(
    rate_per_day: float,
    model_time_unit_days: float,
) -> float:
    """Convert a rate per day to a rate per dimensionless model-time unit.

    This is a unit conversion only. It does not establish that a fitted
    whole-tumour net growth rate equals ModelParameters.cancer_growth, because
    the PDE reaction also contains oxygen/nutrient limitation, immune/drug
    killing and baseline-loss terms.
    """
    if not np.isfinite([rate_per_day, model_time_unit_days]).all():
        raise ValueError("Rates and time scale must be finite.")
    if rate_per_day <= 0 or model_time_unit_days <= 0:
        raise ValueError("Rates and time scale must be positive.")
    return float(rate_per_day * model_time_unit_days)
