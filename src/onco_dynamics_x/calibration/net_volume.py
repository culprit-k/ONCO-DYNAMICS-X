"""Net tumour-volume change estimation for treatment-exposed longitudinal series.

The estimated log-volume slope is a phenomenological net rate. It combines
proliferation, treatment effect, cell death, imaging/segmentation variability
and other mechanisms. It is not an untreated growth rate, drug-killing rate,
EC50, or a calibrated PDE coefficient.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class NetVolumeRateFit:
    """Log-linear net volume change fit for one longitudinal series."""

    net_log_rate_per_day: float
    confidence_interval_95_per_day: list[float]
    standard_error_per_day: float
    intercept_log_volume: float
    r_squared: float
    residual_rmse_log: float
    n_observations: int
    doubling_time_days: float | None
    halving_time_days: float | None
    time_unit: str = "day"
    volume_unit: str = "cc"
    interpretation: str = (
        "treatment-exposed net log-volume slope; not a proliferation or "
        "drug-killing parameter"
    )

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def fit_treatment_exposed_net_volume_rate(
    time_days: NDArray[np.float64] | list[float],
    volume_values: NDArray[np.float64] | list[float],
    *,
    volume_unit: str = "cc",
) -> NetVolumeRateFit:
    """Fit log(volume) against actual elapsed days, allowing shrinkage.

    At least three unique visits are required to estimate residual uncertainty.
    All volumes must be strictly positive. The volume unit is a label because a
    constant unit conversion does not alter the log-volume slope. A zero or
    below-detection FTV is censored and must not be replaced with a pseudocount.
    """
    try:
        from scipy.stats import linregress
        from scipy.stats import t as student_t
    except ImportError as error:
        raise ImportError(
            "Install the science extra: pip install -e '.[science]'."
        ) from error

    time = np.asarray(time_days, dtype=float)
    volume = np.asarray(volume_values, dtype=float)
    if time.ndim != 1 or volume.ndim != 1 or time.shape != volume.shape:
        raise ValueError("time_days and volume_mm3 must be aligned one-dimensional arrays.")
    if len(time) < 3:
        raise ValueError("At least three longitudinal visits are required for uncertainty.")
    if not np.isfinite(time).all() or not np.isfinite(volume).all():
        raise ValueError("Time and volume observations must be finite.")
    if not isinstance(volume_unit, str) or not volume_unit.strip():
        raise ValueError("volume_unit must be a non-empty unit label.")
    if np.any(volume <= 0):
        raise ValueError(
            "Volumes must be positive; zero/below-detection values require a "
            "left-censored observation model, not a pseudocount."
        )

    order = np.argsort(time)
    time = time[order]
    volume = volume[order]
    if np.any(np.diff(time) <= 0):
        raise ValueError("Time points must be unique and strictly increasing.")

    relative_time = time - time[0]
    regression = linregress(relative_time, np.log(volume))
    slope = float(regression.slope)
    stderr = float(regression.stderr)
    degrees_freedom = len(time) - 2
    critical = float(student_t.ppf(0.975, degrees_freedom))
    ci = [float(slope - critical * stderr), float(slope + critical * stderr)]

    fitted_log_volume = float(regression.intercept) + slope * relative_time
    residuals = np.log(volume) - fitted_log_volume
    rmse = float(np.sqrt(np.mean(residuals**2)))
    doubling = float(np.log(2.0) / slope) if slope > 0 else None
    halving = float(np.log(2.0) / abs(slope)) if slope < 0 else None

    return NetVolumeRateFit(
        net_log_rate_per_day=slope,
        confidence_interval_95_per_day=ci,
        standard_error_per_day=stderr,
        intercept_log_volume=float(regression.intercept),
        volume_unit=volume_unit.strip(),
        r_squared=float(regression.rvalue**2),
        residual_rmse_log=rmse,
        n_observations=int(len(time)),
        doubling_time_days=doubling,
        halving_time_days=halving,
    )


def bootstrap_group_net_rates(
    rates_per_day: NDArray[np.float64] | list[float],
    *,
    n_bootstrap: int = 2000,
    seed: int = 42,
) -> dict[str, Any]:
    """Summarise subject-level slopes with a patient-series bootstrap CI."""
    values = np.asarray(rates_per_day, dtype=float)
    if values.ndim != 1 or len(values) < 2:
        raise ValueError("At least two subject-level rates are required.")
    if not np.isfinite(values).all():
        raise ValueError("rates_per_day must contain only finite values.")
    if n_bootstrap < 100:
        raise ValueError("At least 100 bootstrap replicates are required.")

    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(values), size=(n_bootstrap, len(values)))
    medians = np.median(values[draws], axis=1)
    return {
        "subject_series_count": int(len(values)),
        "median_net_log_rate_per_day": float(np.median(values)),
        "mean_net_log_rate_per_day": float(np.mean(values)),
        "median_rate_bootstrap_95ci_per_day": [
            float(np.quantile(medians, 0.025)),
            float(np.quantile(medians, 0.975)),
        ],
        "interpretation": (
            "cohort distribution of subject-level treatment-exposed net rates; "
            "not a causal treatment effect"
        ),
    }
