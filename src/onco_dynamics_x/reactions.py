"""Local reaction terms for the illustrative tumour-microenvironment model."""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.state import TumourState

FloatArray = NDArray[np.float64]


def reaction_terms(
    state: TumourState,
    parameters: ModelParameters,
    infusion_rate: float = 0.0,
) -> dict[str, FloatArray]:
    """Evaluate dimensionless local reactions.

    infusion_rate is an abstract spatially uniform drug-field input, not a dose.
    """
    parameters.validate()
    state.validate()
    if infusion_rate < 0:
        raise ValueError("infusion_rate must be non-negative.")
    p = parameters
    c, immune, fibro = state.cancer, state.immune, state.fibroblast
    oxygen, nutrient, drug = state.oxygen, state.nutrient, state.drug
    evasion, vessels = state.evasion, state.vasculature

    oxygen_factor = oxygen / (p.oxygen_half_saturation + oxygen)
    nutrient_factor = nutrient / (p.nutrient_half_saturation + nutrient)
    growth = p.cancer_growth * c * (1.0 - c / p.cancer_capacity) * oxygen_factor * nutrient_factor
    effective_immune_killing = p.immune_killing / (1.0 + p.evasion_resistance * evasion)
    immune_kill = effective_immune_killing * immune * c
    drug_fraction = drug / (p.drug_ec50 + drug)
    drug_kill = p.drug_killing * drug_fraction * c
    cancer = growth - immune_kill - drug_kill - p.baseline_cancer_loss * c

    recruitment = p.immune_recruitment * c / (p.immune_recruitment_half + c)
    immune_change = (
        p.immune_source + recruitment - p.immune_decay * immune
        - p.fibroblast_immune_suppression * fibro * immune
        - p.evasion_immune_suppression * evasion * immune
    )
    fibroblast = (
        p.fibroblast_activation * c / (p.fibroblast_activation_half + c)
        - p.fibroblast_decay * fibro
    )
    oxygen_change = (
        -p.oxygen_consumption * c * oxygen / (p.oxygen_half_saturation + oxygen)
        + p.oxygen_supply * (p.oxygen_reference * vessels - oxygen)
    )
    nutrient_change = (
        -p.nutrient_consumption * c * nutrient / (p.nutrient_half_saturation + nutrient)
        + p.nutrient_supply * (p.nutrient_reference - nutrient)
    )
    drug_change = infusion_rate - p.drug_decay * drug
    evasion_change = (
        p.evasion_production * c / (p.evasion_production_half + c)
        - p.evasion_decay * evasion
    )
    hypoxia_signal = np.maximum(p.hypoxia_reference - oxygen, 0.0)
    vascular_change = (
        p.vascular_adaptation * hypoxia_signal
        - p.vascular_return * (vessels - 1.0)
    )

    return {
        "cancer": cancer,
        "immune": immune_change,
        "fibroblast": fibroblast,
        "oxygen": oxygen_change,
        "nutrient": nutrient_change,
        "drug": drug_change,
        "evasion": evasion_change,
        "vasculature": vascular_change,
    }
