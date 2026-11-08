"""A small stochastic cell-agent model coupled to PDE fields as an environment."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.state import TumourState

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
CANCER = 0
IMMUNE = 1
FIBROBLAST = 2
_ALLOWED_TYPES = {CANCER, IMMUNE, FIBROBLAST}


@dataclass
class AgentPopulation:
    """Cell positions and phenotypes in the dimensionless square [-1, 1]^2."""

    x: FloatArray
    y: FloatArray
    cell_type: IntArray
    age: FloatArray

    def validate(self) -> None:
        lengths = {len(self.x), len(self.y), len(self.cell_type), len(self.age)}
        if len(lengths) != 1:
            raise ValueError("Agent arrays must have equal lengths.")
        if not (
            np.isfinite(self.x).all()
            and np.isfinite(self.y).all()
            and np.isfinite(self.age).all()
        ):
            raise ValueError("Agent coordinates and ages must be finite.")
        if np.any((self.x < -1.0) | (self.x > 1.0)):
            raise ValueError("Agent x coordinates must lie in [-1, 1].")
        if np.any((self.y < -1.0) | (self.y > 1.0)):
            raise ValueError("Agent y coordinates must lie in [-1, 1].")
        if np.any(self.age < 0):
            raise ValueError("Agent ages must be non-negative.")
        if not set(np.unique(self.cell_type)).issubset(_ALLOWED_TYPES):
            raise ValueError("Unknown cell phenotype.")

    def count(self, cell_type: int) -> int:
        return int(np.count_nonzero(self.cell_type == cell_type))

    def copy(self) -> "AgentPopulation":
        return AgentPopulation(
            self.x.copy(), self.y.copy(), self.cell_type.copy(), self.age.copy()
        )


def initialize_agents(
    cancer_count: int = 80,
    immune_count: int = 50,
    fibroblast_count: int = 35,
    seed: int = 42,
) -> AgentPopulation:
    """Build a deterministic synthetic cell population around a central focus."""
    counts = (cancer_count, immune_count, fibroblast_count)
    if any(count < 0 for count in counts):
        raise ValueError("Agent counts must be non-negative.")
    rng = np.random.default_rng(seed)
    total = sum(counts)
    x = np.empty(total, dtype=float)
    y = np.empty(total, dtype=float)
    cell_type = np.concatenate(
        [
            np.full(cancer_count, CANCER, dtype=np.int64),
            np.full(immune_count, IMMUNE, dtype=np.int64),
            np.full(fibroblast_count, FIBROBLAST, dtype=np.int64),
        ]
    )
    cursor = 0
    for count, radius in (
        (cancer_count, 0.25),
        (immune_count, 0.75),
        (fibroblast_count, 0.55),
    ):
        x[cursor : cursor + count] = rng.normal(0.0, radius / 2.0, count)
        y[cursor : cursor + count] = rng.normal(0.0, radius / 2.0, count)
        cursor += count
    population = AgentPopulation(
        np.clip(x, -1.0, 1.0),
        np.clip(y, -1.0, 1.0),
        cell_type,
        np.zeros(total, dtype=float),
    )
    population.validate()
    return population


def _sample(field: FloatArray, x: FloatArray, y: FloatArray) -> FloatArray:
    """Nearest-grid sampling for positions in the model's [-1, 1] domain."""
    rows = np.rint((y + 1.0) * 0.5 * (field.shape[0] - 1)).astype(int)
    cols = np.rint((x + 1.0) * 0.5 * (field.shape[1] - 1)).astype(int)
    return field[rows, cols]


def _reflect(values: FloatArray) -> FloatArray:
    """Reflect arbitrary overshoots back into the closed unit interval."""
    reflected = np.mod(values + 1.0, 4.0)
    return np.where(reflected > 2.0, 3.0 - reflected, reflected - 1.0)


def _append_agents(
    population: AgentPopulation,
    x: FloatArray,
    y: FloatArray,
    cell_type: int,
    maximum_agents: int,
    rng: np.random.Generator,
) -> AgentPopulation:
    """Append a capped number of new agents at perturbed seed positions."""
    available = max(0, maximum_agents - len(population.x))
    count = min(available, len(x))
    if count == 0:
        return population
    new_x = np.clip(x[:count] + rng.normal(0.0, 0.012, count), -1.0, 1.0)
    new_y = np.clip(y[:count] + rng.normal(0.0, 0.012, count), -1.0, 1.0)
    return AgentPopulation(
        np.concatenate([population.x, new_x]),
        np.concatenate([population.y, new_y]),
        np.concatenate(
            [population.cell_type, np.full(count, cell_type, dtype=np.int64)]
        ),
        np.concatenate([population.age, np.zeros(count)]),
    )


def advance_agents(
    population: AgentPopulation,
    environment: TumourState,
    parameters: ModelParameters,
    dt: float,
    *,
    seed: int = 42,
    maximum_agents: int = 20_000,
) -> AgentPopulation:
    """Advance agents by one stochastic step under local field-dependent rules.

    Cancer agents divide and die, immune agents are recruited and die, and
    fibroblast agents move in a stromal-density-dependent environment. Immune
    chemotaxis is approximated from the cancer-field gradient. This illustrative
    ABM is not a validated model of cell kinetics.
    """
    if dt <= 0:
        raise ValueError("dt must be positive.")
    if maximum_agents < len(population.x):
        raise ValueError("maximum_agents cannot be lower than the current population.")
    population.validate()
    environment.validate()
    parameters.validate()
    rng = np.random.default_rng(seed)
    result = population.copy()
    result.age += dt

    cancer_mask = result.cell_type == CANCER
    immune_mask = result.cell_type == IMMUNE
    fibro_mask = result.cell_type == FIBROBLAST

    cancer_field = _sample(environment.cancer, result.x, result.y)
    oxygen = _sample(environment.oxygen, result.x, result.y)
    nutrient = _sample(environment.nutrient, result.x, result.y)
    local_immune = _sample(environment.immune, result.x, result.y)
    drug = _sample(environment.drug, result.x, result.y)
    fibro_field = _sample(environment.fibroblast, result.x, result.y)
    evasion = _sample(environment.evasion, result.x, result.y)

    motility = np.where(
        cancer_mask, 0.008, np.where(immune_mask, 0.025, 0.006 / (1.0 + fibro_field))
    )
    grad_y, grad_x = np.gradient(
        environment.cancer, 2.0 / (environment.cancer.shape[0] - 1)
    )
    immune_chemotaxis_x = _sample(grad_x, result.x, result.y)
    immune_chemotaxis_y = _sample(grad_y, result.x, result.y)
    result.x += (
        rng.normal(0.0, 1.0, len(result.x)) * np.sqrt(2.0 * motility * dt)
        + np.where(immune_mask, 0.01 * immune_chemotaxis_x * dt, 0.0)
    )
    result.y += (
        rng.normal(0.0, 1.0, len(result.y)) * np.sqrt(2.0 * motility * dt)
        + np.where(immune_mask, 0.01 * immune_chemotaxis_y * dt, 0.0)
    )
    result.x = _reflect(result.x)
    result.y = _reflect(result.y)

    oxygen_factor = oxygen / (parameters.oxygen_half_saturation + oxygen)
    nutrient_factor = nutrient / (parameters.nutrient_half_saturation + nutrient)
    immune_hazard = (
        parameters.immune_killing
        * local_immune
        / (1.0 + parameters.evasion_resistance * evasion)
    )
    drug_hazard = parameters.drug_killing * drug / (parameters.drug_ec50 + drug)
    cancer_hazard = parameters.baseline_cancer_loss + immune_hazard + drug_hazard
    cancer_birth_hazard = parameters.cancer_growth * oxygen_factor * nutrient_factor

    death_probability = 1.0 - np.exp(-cancer_hazard * dt)
    birth_probability = 1.0 - np.exp(-cancer_birth_hazard * dt)
    deaths = cancer_mask & (rng.random(len(result.x)) < death_probability)
    daughter_mask = (
        cancer_mask
        & ~deaths
        & (rng.random(len(result.x)) < birth_probability)
    )
    daughter_indices = np.flatnonzero(daughter_mask)

    immune_recruitment = (
        parameters.immune_recruitment
        * cancer_field
        / (parameters.immune_recruitment_half + cancer_field)
    )
    recruit_immune_mask = (
        rng.random(len(result.x)) < (1.0 - np.exp(-immune_recruitment * dt))
    )
    fibroblast_recruitment = (
        parameters.fibroblast_activation
        * cancer_field
        / (parameters.fibroblast_activation_half + cancer_field)
    )
    recruit_fibro_mask = (
        rng.random(len(result.x)) < (1.0 - np.exp(-fibroblast_recruitment * dt))
    )

    immune_death_probability = 1.0 - np.exp(-parameters.immune_decay * dt)
    fibro_death_probability = 1.0 - np.exp(-parameters.fibroblast_decay * dt)
    immune_deaths = immune_mask & (
        rng.random(len(result.x)) < immune_death_probability
    )
    fibro_deaths = fibro_mask & (
        rng.random(len(result.x)) < fibro_death_probability
    )

    daughter_x, daughter_y = result.x[daughter_indices], result.y[daughter_indices]
    immune_x = result.x[recruit_immune_mask]
    immune_y = result.y[recruit_immune_mask]
    fibro_x = result.x[recruit_fibro_mask]
    fibro_y = result.y[recruit_fibro_mask]

    keep = ~(deaths | immune_deaths | fibro_deaths)
    result = AgentPopulation(
        x=result.x[keep],
        y=result.y[keep],
        cell_type=result.cell_type[keep],
        age=result.age[keep],
    )
    result = _append_agents(
        result, daughter_x, daughter_y, CANCER, maximum_agents, rng
    )
    result = _append_agents(
        result, immune_x, immune_y, IMMUNE, maximum_agents, rng
    )
    result = _append_agents(
        result, fibro_x, fibro_y, FIBROBLAST, maximum_agents, rng
    )
    result.validate()
    return result
