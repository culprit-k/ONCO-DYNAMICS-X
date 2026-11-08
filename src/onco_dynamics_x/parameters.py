"""Dimensionless parameters for the illustrative reaction-diffusion model."""
from __future__ import annotations

from dataclasses import dataclass, fields


@dataclass(frozen=True)
class ModelParameters:
    """Provisional dimensionless parameters; none are clinically calibrated."""

    diffusion_cancer: float = 0.010
    diffusion_immune: float = 0.018
    diffusion_fibroblast: float = 0.006
    diffusion_oxygen: float = 0.040
    diffusion_nutrient: float = 0.030
    diffusion_drug: float = 0.035
    diffusion_evasion: float = 0.008
    diffusion_vasculature: float = 0.002

    cancer_growth: float = 0.55
    cancer_capacity: float = 1.0
    immune_killing: float = 0.30
    baseline_cancer_loss: float = 0.025
    drug_killing: float = 0.50
    drug_ec50: float = 0.25
    evasion_resistance: float = 1.5
    oxygen_half_saturation: float = 0.20
    nutrient_half_saturation: float = 0.20

    immune_source: float = 0.015
    immune_recruitment: float = 0.10
    immune_recruitment_half: float = 0.20
    immune_decay: float = 0.12
    fibroblast_immune_suppression: float = 0.16
    evasion_immune_suppression: float = 0.22

    fibroblast_activation: float = 0.12
    fibroblast_activation_half: float = 0.20
    fibroblast_decay: float = 0.08

    oxygen_consumption: float = 0.16
    oxygen_supply: float = 0.30
    oxygen_reference: float = 1.0
    nutrient_consumption: float = 0.11
    nutrient_supply: float = 0.22
    nutrient_reference: float = 1.0

    drug_decay: float = 0.18
    evasion_production: float = 0.08
    evasion_production_half: float = 0.20
    evasion_decay: float = 0.09

    hypoxia_reference: float = 0.65
    vascular_adaptation: float = 0.025
    vascular_return: float = 0.035

    maximum_drug: float = 5.0
    maximum_field: float = 2.0

    def validate(self) -> None:
        for parameter in fields(self):
            if getattr(self, parameter.name) < 0:
                raise ValueError(f"{parameter.name} must be non-negative.")
        if self.cancer_capacity <= 0 or self.drug_ec50 <= 0:
            raise ValueError("Cancer capacity and drug_ec50 must be positive.")
        if self.oxygen_half_saturation <= 0 or self.nutrient_half_saturation <= 0:
            raise ValueError("Oxygen and nutrient half-saturation values must be positive.")
        if self.maximum_drug <= 0 or self.maximum_field <= 0:
            raise ValueError("Numerical clipping limits must be positive.")

    @property
    def diffusion_coefficients(self) -> dict[str, float]:
        return {
            "cancer": self.diffusion_cancer,
            "immune": self.diffusion_immune,
            "fibroblast": self.diffusion_fibroblast,
            "oxygen": self.diffusion_oxygen,
            "nutrient": self.diffusion_nutrient,
            "drug": self.diffusion_drug,
            "evasion": self.diffusion_evasion,
            "vasculature": self.diffusion_vasculature,
        }
