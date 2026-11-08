"""Validated local serial protocol for the educational hardware console."""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, replace
from typing import Any

from onco_dynamics_x.parameters import ModelParameters

_ALLOWED_KEYS = {
    "event",
    "tumour_growth",
    "immune_activity",
    "run",
    "reset",
}


@dataclass(frozen=True)
class ControlMessage:
    """Validated controls normalized to [0, 1]; not biological measurements."""

    tumour_growth: float
    immune_activity: float
    run: bool = False
    reset: bool = False


def _normalized_value(payload: dict[str, Any], key: str) -> float:
    value = payload.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"'{key}' must be a numeric value.")
    parsed = float(value)
    if not math.isfinite(parsed) or not 0.0 <= parsed <= 1.0:
        raise ValueError(f"'{key}' must be finite and within [0, 1].")
    return parsed


def parse_control_message(line: str) -> ControlMessage:
    """Parse one JSON line and reject unknown keys or invalid control ranges."""
    try:
        payload = json.loads(line)
    except json.JSONDecodeError as error:
        raise ValueError("Serial message is not valid JSON.") from error
    if not isinstance(payload, dict):
        raise ValueError("Serial message must be a JSON object.")
    unknown = set(payload) - _ALLOWED_KEYS
    if unknown:
        raise ValueError(f"Unknown message key(s): {', '.join(sorted(unknown))}.")
    if payload.get("event") != "controls":
        raise ValueError("Expected event='controls'.")
    for key in ("run", "reset"):
        if key in payload and not isinstance(payload[key], bool):
            raise ValueError(f"'{key}' must be a boolean.")
    return ControlMessage(
        tumour_growth=_normalized_value(payload, "tumour_growth"),
        immune_activity=_normalized_value(payload, "immune_activity"),
        run=payload.get("run", False),
        reset=payload.get("reset", False),
    )


def apply_controls(
    parameters: ModelParameters,
    message: ControlMessage,
) -> ModelParameters:
    """Map normalized knobs to bounded model parameters for demonstration use."""
    updated = replace(
        parameters,
        cancer_growth=0.10 + 1.10 * message.tumour_growth,
        immune_killing=0.05 + 0.75 * message.immune_activity,
    )
    updated.validate()
    return updated
