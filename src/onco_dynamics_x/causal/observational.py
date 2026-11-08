"""Optional observational causal-analysis wrapper with explicit assumptions."""
from __future__ import annotations

from typing import Any

import pandas as pd


def estimate_backdoor_linear_effect(
    data: pd.DataFrame,
    treatment: str,
    outcome: str,
    common_causes: list[str],
) -> tuple[Any, Any]:
    """Identify and estimate a linear backdoor effect using declared confounders.

    This wrapper does not discover a valid causal graph. The analyst must justify
    treatment, outcome, common causes, temporal ordering, selection assumptions,
    and positivity before interpreting an estimate causally.
    """
    try:
        from dowhy import CausalModel
    except ImportError as error:
        raise ImportError(
            "Install the causal extra with pip install -e .[causal]."
        ) from error

    required = [treatment, outcome, *common_causes]
    if not data.index.is_unique:
        raise ValueError("DataFrame index must be unique.")
    missing = [column for column in required if column not in data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}.")
    if not common_causes:
        raise ValueError("Declare at least one common cause; do not assume none.")
    if data[required].isna().any().any():
        raise ValueError("Required causal-analysis columns contain missing values.")
    if data.empty:
        raise ValueError("Causal-analysis data cannot be empty.")

    model = CausalModel(
        data=data,
        treatment=treatment,
        outcome=outcome,
        common_causes=common_causes,
    )
    estimand = model.identify_effect(proceed_when_unidentifiable=False)
    estimate = model.estimate_effect(
        estimand,
        method_name="backdoor.linear_regression",
        confidence_intervals=True,
    )
    return estimand, estimate
