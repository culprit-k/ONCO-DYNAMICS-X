"""Optional PyMC parameter posterior construction for user-supplied forward models."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import numpy as np


@dataclass(frozen=True)
class PriorSpec:
    """Scalar prior specification for an unknown model parameter."""

    distribution: str
    location: float
    scale: float

    def validate(self) -> None:
        if self.distribution not in {"normal", "lognormal"}:
            raise ValueError("Prior distribution must be 'normal' or 'lognormal'.")
        if self.scale <= 0 or not np.isfinite(self.scale):
            raise ValueError("Prior scale must be finite and positive.")
        if not np.isfinite(self.location):
            raise ValueError("Prior location must be finite.")


def build_gaussian_parameter_model(
    observations: np.ndarray,
    forward_model: Callable[[dict[str, Any]], Any],
    priors: dict[str, PriorSpec],
    observation_sd: float,
) -> Any:
    """Build a Gaussian likelihood model around a backend-compatible predictor.

    forward_model must operate on PyTensor/PyMC random variables and return an
    expected-observation tensor compatible with observations. Scientific priors
    must be selected by the analyst; this utility does not infer biological priors.
    """
    observed = np.asarray(observations, dtype=float)
    if observed.size == 0 or not np.isfinite(observed).all():
        raise ValueError("observations must be non-empty and finite.")
    if observation_sd <= 0 or not np.isfinite(observation_sd):
        raise ValueError("observation_sd must be finite and positive.")
    if not priors:
        raise ValueError("At least one parameter prior is required.")
    for name, prior in priors.items():
        if not name.isidentifier():
            raise ValueError(f"Invalid parameter name: {name}")
        prior.validate()
    try:
        import pymc as pm
    except ImportError as error:
        raise ImportError(
            "Install the Bayesian extra with pip install -e .[bayesian]."
        ) from error

    with pm.Model() as model:
        parameters: dict[str, Any] = {}
        for name, prior in priors.items():
            if prior.distribution == "normal":
                parameters[name] = pm.Normal(
                    name, mu=prior.location, sigma=prior.scale
                )
            else:
                parameters[name] = pm.LogNormal(
                    name, mu=prior.location, sigma=prior.scale
                )
        expected = forward_model(parameters)
        pm.Normal(
            "observations",
            mu=expected,
            sigma=observation_sd,
            observed=observed,
        )
    return model


def sample_parameter_posterior(
    model: Any,
    draws: int = 1_000,
    tune: int = 1_000,
    chains: int = 4,
    random_seed: int = 42,
) -> Any:
    """Sample a parameter posterior from a model built by this module."""
    if min(draws, tune, chains) < 1:
        raise ValueError("draws, tune and chains must all be positive.")
    try:
        import pymc as pm
    except ImportError as error:
        raise ImportError(
            "Install the Bayesian extra with pip install -e .[bayesian]."
        ) from error
    with model:
        return pm.sample(
            draws=draws,
            tune=tune,
            chains=chains,
            random_seed=random_seed,
            return_inferencedata=True,
        )
