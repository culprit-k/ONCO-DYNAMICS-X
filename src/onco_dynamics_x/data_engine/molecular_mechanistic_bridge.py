"""Exploratory bridge from fixed-panel expression to model-facing descriptors.

The scores in this module are descriptive summaries of rank-normalised
expression. They are not pathway-activity estimates, calibrated rates, or
patient-specific treatment-response predictions. No PDE/ABM parameter is
modified by this module.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]

# Prespecified, interpretable groupings drawn from the fixed research panel.
# These are heuristic feature modules, not validated gene-set signatures.
MOLECULAR_MODULES: dict[str, tuple[str, ...]] = {
    "proliferation_expression": (
        "ANLN", "BIRC5", "CCNB1", "CDC20", "CENPF", "CEP55",
        "KIF2C", "MKI67", "MYBL2", "TYMS", "UBE2C",
    ),
    "estrogen_luminal_expression": (
        "BAG1", "BCL2", "ESR1", "FOXA1", "MAPT", "MLPH",
        "NAT1", "PGR", "SLC39A6",
    ),
    "her2_associated_expression": ("ERBB2", "GRB7"),
    "basal_epithelial_expression": (
        "CDH3", "EGFR", "FOXC1", "KRT14", "KRT17", "KRT5",
        "MMP11", "TMEM45B",
    ),
}

# These are hypotheses for future calibration, not executable parameter maps.
MECHANISTIC_LINK_HYPOTHESES: dict[str, dict[str, Any]] = {
    "proliferation_expression": {
        "candidate_model_target": "cancer-cell proliferation propensity",
        "possible_model_component": "ABM phenotype transition or proliferation-rate prior",
        "evidence_needed_before_use": [
            "independent proliferation measurement such as Ki-67 or validated cell-cycle assay",
            "patient- or specimen-matched growth/kinetic observations",
            "external calibration and uncertainty analysis",
        ],
    },
    "estrogen_luminal_expression": {
        "candidate_model_target": "luminal/estrogen-associated phenotype state",
        "possible_model_component": "cell-state composition or phenotype-transition prior",
        "evidence_needed_before_use": [
            "ER/PR immunohistochemistry or receptor status with provenance",
            "independent cohort replication",
            "a prespecified mapping model with uncertainty",
        ],
    },
    "her2_associated_expression": {
        "candidate_model_target": "HER2-associated molecular state",
        "possible_model_component": "phenotype-state annotation only",
        "evidence_needed_before_use": [
            "HER2 IHC/ISH or copy-number amplification evidence",
            "ideally receptor/phosphoprotein measurements",
            "external calibration; expression alone is not receptor activity",
        ],
    },
    "basal_epithelial_expression": {
        "candidate_model_target": "basal/epithelial phenotype state",
        "possible_model_component": "cell-state composition or phenotype-transition prior",
        "evidence_needed_before_use": [
            "independent basal-marker pathology or single-cell/spatial evidence",
            "external cohort replication",
            "a prespecified mapping model with uncertainty",
        ],
    },
}


@dataclass(frozen=True)
class MolecularMechanisticProfile:
    """Sample-level descriptive module scores; no mechanistic parameters."""

    sample_scores: FloatArray
    module_names: tuple[str, ...]
    module_genes: dict[str, tuple[str, ...]]
    gene_coverage: dict[str, float]
    interpretation: str = (
        "mean within-cohort rank-normal scores; descriptive expression summaries, "
        "not pathway activity or calibrated biological rates"
    )


def score_molecular_modules(
    expression_scores: NDArray[np.float64],
    gene_symbols: tuple[str, ...] | list[str],
    *,
    minimum_gene_coverage: float = 0.75,
) -> MolecularMechanisticProfile:
    """Compute outcome-blind module means from rank-normalised expression.

    Input columns must already be aligned to unique gene symbols and transformed
    without outcome/subtype labels. Module means are not re-standardised, and
    the function does not fit thresholds or alter any simulation parameter.
    """
    matrix = np.asarray(expression_scores, dtype=float)
    genes = tuple(str(gene).upper() for gene in gene_symbols)
    if matrix.ndim != 2 or matrix.shape[1] != len(genes):
        raise ValueError("expression_scores must be [samples, genes] and match gene_symbols.")
    if matrix.shape[0] < 1:
        raise ValueError("At least one sample is required.")
    if not np.isfinite(matrix).all():
        raise ValueError("expression_scores must contain only finite values.")
    if len(set(genes)) != len(genes):
        raise ValueError("gene_symbols must be unique.")
    if not 0 < minimum_gene_coverage <= 1:
        raise ValueError("minimum_gene_coverage must lie in (0, 1].")

    gene_to_index = {gene: index for index, gene in enumerate(genes)}
    scores: list[FloatArray] = []
    coverage: dict[str, float] = {}
    used_genes: dict[str, tuple[str, ...]] = {}
    for module_name, module_genes in MOLECULAR_MODULES.items():
        present = tuple(gene for gene in module_genes if gene in gene_to_index)
        fraction = len(present) / len(module_genes)
        coverage[module_name] = float(fraction)
        if fraction < minimum_gene_coverage:
            raise ValueError(
                f"{module_name} gene coverage {fraction:.1%} is below "
                f"{minimum_gene_coverage:.1%}."
            )
        indices = [gene_to_index[gene] for gene in present]
        scores.append(np.mean(matrix[:, indices], axis=1))
        used_genes[module_name] = present

    return MolecularMechanisticProfile(
        sample_scores=np.column_stack(scores),
        module_names=tuple(MOLECULAR_MODULES),
        module_genes=used_genes,
        gene_coverage=coverage,
    )


def describe_mechanistic_link_policy() -> dict[str, Any]:
    """Return explicit, non-executable hypotheses and forbidden shortcuts."""
    return {
        "status": "hypothesis registry only; no parameter updates",
        "module_hypotheses": MECHANISTIC_LINK_HYPOTHESES,
        "prohibited_direct_mappings": [
            "expression score directly to PDE diffusion coefficient",
            "expression score directly to oxygen consumption or hypoxia rate",
            "expression score directly to drug sensitivity or treatment response",
            "expression score directly to immune killing or immune-evasion rate",
            "expression score directly to a patient-specific tumour growth rate",
        ],
        "required_calibration_sequence": [
            "define a measurable biological target for each proposed parameter",
            "obtain matched molecular and biological measurements",
            "fit mapping only in a designated training cohort",
            "quantify identifiability, uncertainty and sensitivity",
            "freeze mapping and evaluate once in an independent cohort",
            "retain a no-molecular-input baseline and report failure cases",
        ],
    }
