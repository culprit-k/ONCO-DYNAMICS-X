"""Extended-Tofts and spoiled-GRE MRI observation components.

This module does not prescribe an unvalidated mapping from PDE cell density to
contrast kinetics. A future learned field-to-parameter network must produce
constrained Ktrans, ve and vp maps, which are then passed through this forward
MRI signal model using acquisition-specific inputs.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def constrain_pk_parameter_logits(logits: np.ndarray) -> dict[str, np.ndarray]:
    """Map [..., 3] unconstrained logits to physical PK parameter ranges.

    The last axis is ordered as Ktrans (min^-1), ve and vp. Ktrans is softplus
    transformed; ve is constrained to (0,1); vp is constrained to (0,1-ve).
    The neural network producing these logits must be learned on training
    patients only. This function contains no learned coefficients.
    """
    values = np.asarray(logits, dtype=float)
    if values.ndim < 1 or values.shape[-1] != 3:
        raise ValueError("logits must have a final dimension of length 3")
    if not np.isfinite(values).all():
        raise ValueError("logits must be finite")
    clipped = np.clip(values, -30.0, 30.0)
    ktrans = np.logaddexp(0.0, clipped[..., 0])
    ve = 1.0 / (1.0 + np.exp(-clipped[..., 1]))
    vp_fraction = 1.0 / (1.0 + np.exp(-clipped[..., 2]))
    vp = (1.0 - ve) * vp_fraction
    return {
        "ktrans_per_min": ktrans,
        "ve_fraction": ve,
        "vp_fraction": vp,
    }


def extended_tofts_concentration(
    times_seconds: np.ndarray,
    arterial_input_mM: np.ndarray,
    *,
    ktrans_per_min: np.ndarray | float,
    ve_fraction: np.ndarray | float,
    vp_fraction: np.ndarray | float,
) -> np.ndarray:
    """Calculate tissue contrast concentration with the extended Tofts model.

    The time grid is actual elapsed time after contrast injection, expressed in
    seconds; Ktrans is min^-1; concentrations are mmol/L (mM). The numerical
    convolution uses a recursive trapezoidal approximation. Parameter maps are
    broadcast over spatial voxels and time is the leading output axis.
    """
    time_s = np.asarray(times_seconds, dtype=float)
    cp = np.asarray(arterial_input_mM, dtype=float)
    if time_s.ndim != 1 or cp.ndim != 1 or time_s.shape != cp.shape or time_s.size < 2:
        raise ValueError("times_seconds and arterial_input_mM must be equal-length 1D arrays")
    if not np.isfinite(time_s).all() or not np.isfinite(cp).all():
        raise ValueError("Time and arterial input arrays must be finite")
    if time_s[0] < 0 or np.any(np.diff(time_s) <= 0):
        raise ValueError("times_seconds must start at >=0 and be strictly increasing")
    if (cp < 0).any():
        raise ValueError("Arterial input concentration cannot be negative")

    ktrans, ve, vp = np.broadcast_arrays(
        np.asarray(ktrans_per_min, dtype=float),
        np.asarray(ve_fraction, dtype=float),
        np.asarray(vp_fraction, dtype=float),
    )
    if not np.isfinite(ktrans).all() or not np.isfinite(ve).all() or not np.isfinite(vp).all():
        raise ValueError("PK parameters must be finite")
    if (ktrans < 0).any():
        raise ValueError("Ktrans must be non-negative")
    if ((ve <= 0) | (ve >= 1)).any():
        raise ValueError("ve must be strictly between 0 and 1")
    if ((vp < 0) | (vp >= 1) | (vp + ve > 1.0 + 1e-10)).any():
        raise ValueError("vp must be non-negative and vp + ve cannot exceed 1")

    kep = ktrans / ve
    integral = np.zeros_like(kep, dtype=float)
    result = np.empty((time_s.size,) + ktrans.shape, dtype=float)
    result[0] = vp * cp[0]
    for index in range(1, time_s.size):
        dt_min = (time_s[index] - time_s[index - 1]) / 60.0
        decay = np.exp(-kep * dt_min)
        integral = (
            decay * integral
            + 0.5 * dt_min * (cp[index - 1] * decay + cp[index])
        )
        result[index] = vp * cp[index] + ktrans * integral
    return result


def spoiled_gre_percent_enhancement(
    tissue_concentration_mM: np.ndarray,
    *,
    baseline_t1_seconds: np.ndarray | float,
    flip_angle_degrees: float,
    repetition_time_seconds: float,
    longitudinal_relaxivity_s_inv_per_mM: float,
) -> np.ndarray:
    """Convert concentration to SPGR percent enhancement under a T2*-free model.

    Uses the spoiled-GRE steady-state signal equation with actual baseline T1,
    flip angle, TR and relaxivity. Coil gain/M0 cancels in the signal ratio.
    T2* attenuation, B1 inhomogeneity, inflow and motion are not represented;
    these must be corrected/modelled or included in an uncertainty analysis
    before interpreting real patient maps.
    """
    concentration = np.asarray(tissue_concentration_mM, dtype=float)
    if concentration.ndim < 1 or not np.isfinite(concentration).all():
        raise ValueError("tissue_concentration_mM must be a finite array with time first")
    if (concentration < 0).any():
        raise ValueError("Tissue contrast concentration cannot be negative")
    t10 = np.asarray(baseline_t1_seconds, dtype=float)
    if not np.isfinite(t10).all() or (t10 <= 0).any():
        raise ValueError("baseline_t1_seconds must be finite and positive")
    alpha_deg = float(flip_angle_degrees)
    tr = float(repetition_time_seconds)
    relaxivity = float(longitudinal_relaxivity_s_inv_per_mM)
    if not np.isfinite(alpha_deg) or not 0 < alpha_deg < 180:
        raise ValueError("flip_angle_degrees must be between 0 and 180")
    if not np.isfinite(tr) or tr <= 0:
        raise ValueError("repetition_time_seconds must be finite and positive")
    if not np.isfinite(relaxivity) or relaxivity < 0:
        raise ValueError("longitudinal_relaxivity must be finite and non-negative")

    try:
        np.broadcast_shapes(concentration.shape[1:], t10.shape)
    except ValueError as exc:
        raise ValueError("Baseline T1 map is not broadcast-compatible with concentration voxels") from exc

    alpha = np.deg2rad(alpha_deg)
    r10 = 1.0 / t10
    r1 = r10 + relaxivity * concentration
    e1 = np.exp(-tr * r1)
    e10 = np.exp(-tr * r10)
    signal = np.sin(alpha) * (1.0 - e1) / (1.0 - np.cos(alpha) * e1)
    baseline_signal = np.sin(alpha) * (1.0 - e10) / (
        1.0 - np.cos(alpha) * e10
    )
    enhancement = 100.0 * (signal / baseline_signal - 1.0)
    return enhancement
