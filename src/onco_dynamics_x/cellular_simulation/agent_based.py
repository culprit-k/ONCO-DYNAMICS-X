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

    This illustrative ABM is not a validated model of cell kinetics.
    Daughter cells are placed near parents. Boundaries are reflecting.
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

    oxygen = _sample(environment.oxygen, result.x, result.y)
    nutrient = _sample(environment.nutrient, result.x, result.y)
    local_immune = _sample(environment.immune, result.x, result.y)
    drug = _sample(environment.drug, result.x, result.y)
    evasion = _sample(environment.evasion, result.x, result.y)

    motility = np.where(cancer_mask, 0.008, np.where(immune_mask, 0.025, 0.006))
    result.x += rng.normal(0.0, 1.0, len(result.x)) * np.sqrt(2.0 * motility * dt)
    result.y += rng.normal(0.0, 1.0, len(result.y)) * np.sqrt(2.0 * motility * dt)
    result.x = np.where(result.x > 1.0, 2.0 - result.x, result.x)
    result.x = np.where(result.x < -1.0, -2.0 - result.x, result.x)
    result.y = np.where(result.y > 1.0, 2.0 - result.y, result.y)
    result.y = np.where(result.y < -1.0, -2.0 - result.y, result.y)

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
    random_values = rng.random(len(result.x))
    deaths = cancer_mask & (random_values < death_probability)
    daughter_mask = (
        cancer_mask
        & ~deaths
        & (rng.random(len(result.x)) < birth_probability)
    )
    daughter_indices = np.flatnonzero(daughter_mask)
    if len(result.x) + len(daughter_indices) > maximum_agents:
        allowed = max(0, maximum_agents - len(result.x))
        daughter_indices = daughter_indices[:allowed]

    immune_death_probability = 1.0 - np.exp(-parameters.immune_decay * dt)
    fibro_death_probability = 1.0 - np.exp(-parameters.fibroblast_decay * dt)
    immune_deaths = immune_mask & (
        rng.random(len(result.x)) < immune_death_probability
    )
    fibro_deaths = fibro_mask & (
        rng.random(len(result.x)) < fibro_death_probability
    )
    keep = ~(deaths | immune_deaths | fibro_deaths)
    result = AgentPopulation(
        x=result.x[keep],
        y=result.y[keep],
        cell_type=result.cell_type[keep],
        age=result.age[keep],
    )
    if len(daughter_indices):
        daughters_x = (
            population.x[daughter_indices]
            + rng.normal(0.0, 0.015, len(daughter_indices))
        )
        daughters_y = (
            population.y[daughter_indices]
            + rng.normal(0.0, 0.015, len(daughter_indices))
        )
        result.x = np.concatenate([result.x, np.clip(daughters_x, -1.0, 1.0)])
        result.y = np.concatenate([result.y, np.clip(daughters_y, -1.0, 1.0)])
        result.cell_type = np.concatenate(
            [result.cell_type, np.full(len(daughter_indices), CANCER, dtype=np.int64)]
        )
        result.age = np.concatenate([result.age, np.zeros(len(daughter_indices))])

    result.validate()
    return result
