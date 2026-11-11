"""Leakage-conscious research features from the published PAM50 gene panel.

This module deliberately does not implement or reproduce the commercial
Prosigna assay or the proprietary clinical Risk-of-Recurrence score. It uses
the 50-gene research panel as a fixed feature inventory, then applies an
explicit cohort-wise rank-normal transformation for exploratory cross-platform
analysis.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.special import ndtri
from scipy.stats import rankdata
from sklearn.decomposition import PCA

FloatArray = NDArray[np.float64]

# The 50-gene panel as exposed by the genefu PAM50 gene mapping.
# These are gene identifiers only; no subtype centroids or proprietary
# Prosigna coefficients are included in this repository.
PAM50_GENES: tuple[str, ...] = (
    "ACTR3B", "ANLN", "BAG1", "BCL2", "BIRC5", "BLVRA", "CCNB1", "CCNE1",
    "CDC20", "CDC6", "CDCA1", "CDH3", "CENPF", "CEP55", "CXXC5", "EGFR",
    "ERBB2", "ESR1", "EXO1", "FGFR4", "FOXA1", "FOXC1", "GPR160", "GRB7",
    "KIF2C", "KNTC2", "KRT14", "KRT17", "KRT5", "MAPT", "MDM2", "MELK",
    "MIA", "MKI67", "MLPH", "MMP11", "MYBL2", "MYC", "NAT1", "ORC6L",
    "PGR", "PHGDH", "PTTG1", "RRM2", "SFRP1", "SLC39A6", "TMEM45B", "TYMS",
    "UBE2C", "UBE2T",
)

# Historical identifiers retained by the published PAM50 mapping; cBioPortal
# can expose the current HUGO symbol for the same Entrez Gene ID.
PAM50_SYMBOL_ALIASES: dict[str, tuple[str, ...]] = {
    "CDCA1": ("CDCA1", "NUF2"),
    "KNTC2": ("KNTC2", "NDC80"),
    "ORC6L": ("ORC6L", "ORC6"),
}


@dataclass(frozen=True)
class RankNormalisedExpression:
    """Outcome-blind, per-cohort rank-normalised expression matrix."""

    values: FloatArray
    gene_symbols: tuple[str, ...]
    missing_fraction_by_gene: dict[str, float]
    samples_with_imputation: int


def rank_normalise_expression(
    values: NDArray[np.float64],
    gene_symbols: tuple[str, ...] | list[str],
    *,
    minimum_observed_fraction: float = 0.95,
) -> RankNormalisedExpression:
    """Rank-normalise each gene across samples, without using clinical labels.

    For each gene, finite values are replaced by inverse-normal scores of their
    within-cohort mid-ranks. Missing values are set to zero, the cohort median
    on the normal-score scale. This is a cohort-level, transductive transform:
    all expression values in the evaluation cohort are used, but no outcomes,
    subtype labels or model residuals are used.
    """
    matrix = np.asarray(values, dtype=float)
    genes = tuple(str(gene) for gene in gene_symbols)
    if matrix.ndim != 2 or matrix.shape[1] != len(genes):
        raise ValueError("values must be [samples, genes] and match gene_symbols.")
    if matrix.shape[0] < 3:
        raise ValueError("At least three samples are required for rank normalisation.")
    if len(set(genes)) != len(genes):
        raise ValueError("gene_symbols must be unique.")
    if not 0 < minimum_observed_fraction <= 1:
        raise ValueError("minimum_observed_fraction must lie in (0, 1].")

    transformed = np.zeros(matrix.shape, dtype=float)
    missing_fraction: dict[str, float] = {}
    any_imputed = np.zeros(matrix.shape[0], dtype=bool)
    for column, gene in enumerate(genes):
        finite = np.isfinite(matrix[:, column])
        observed_count = int(finite.sum())
        coverage = observed_count / matrix.shape[0]
        missing_fraction[gene] = float(1.0 - coverage)
        if coverage < minimum_observed_fraction:
            raise ValueError(
                f"{gene} observed in only {coverage:.1%} of samples; "
                f"minimum is {minimum_observed_fraction:.1%}."
            )
        if observed_count < 3:
            raise ValueError(f"{gene} has fewer than three finite measurements.")
        ranks = rankdata(matrix[finite, column], method="average")
        probabilities = (ranks - 0.5) / observed_count
        transformed[finite, column] = ndtri(
            np.clip(probabilities, 1e-6, 1.0 - 1e-6)
        )
        any_imputed |= ~finite
    return RankNormalisedExpression(
        values=transformed,
        gene_symbols=genes,
        missing_fraction_by_gene=missing_fraction,
        samples_with_imputation=int(any_imputed.sum()),
    )


def fit_tcga_pca_projection(
    training_features: NDArray[np.float64],
    external_features: NDArray[np.float64],
    gene_symbols: tuple[str, ...] | list[str],
    *,
    n_components: int = 5,
    random_state: int = 42,
) -> tuple[FloatArray, FloatArray, dict[str, Any]]:
    """Fit PCA in TCGA and apply the frozen loadings to external features."""
    train = np.asarray(training_features, dtype=float)
    external = np.asarray(external_features, dtype=float)
    genes = tuple(gene_symbols)
    if train.ndim != 2 or external.ndim != 2:
        raise ValueError("Training and external features must be 2D matrices.")
    if train.shape[1] != len(genes) or external.shape[1] != len(genes):
        raise ValueError("Both matrices must use the same ordered gene list.")
    if not np.isfinite(train).all() or not np.isfinite(external).all():
        raise ValueError("PCA input matrices must contain only finite values.")
    if n_components < 1 or n_components > min(train.shape):
        raise ValueError("n_components must be within the training matrix dimensions.")

    model = PCA(n_components=n_components, svd_solver="full", random_state=random_state)
    raw_train_scores = model.fit_transform(train)
    raw_external_scores = model.transform(external)
    score_means = np.mean(raw_train_scores, axis=0)
    score_scales = np.std(raw_train_scores, axis=0, ddof=0)
    if not np.isfinite(score_scales).all() or np.any(score_scales <= 0):
        raise ValueError("A fitted PCA component has zero or invalid training variance.")
    train_scores = (raw_train_scores - score_means) / score_scales
    external_scores = (raw_external_scores - score_means) / score_scales
    loadings = []
    for component_index, vector in enumerate(model.components_):
        positive = np.argsort(vector)[::-1][:6]
        negative = np.argsort(vector)[:6]
        loadings.append(
            {
                "component": f"PC{component_index + 1}",
                "explained_variance_ratio": float(
                    model.explained_variance_ratio_[component_index]
                ),
                "top_positive_loadings": [
                    {"gene": genes[int(index)], "loading": float(vector[index])}
                    for index in positive
                ],
                "top_negative_loadings": [
                    {"gene": genes[int(index)], "loading": float(vector[index])}
                    for index in negative
                ],
            }
        )
    metadata = {
        "fit_cohort": "TCGA-BRCA",
        "external_transform": "frozen TCGA PCA components",
        "n_components": int(n_components),
        "random_state": int(random_state),
        "score_standardisation": "zero mean and unit standard deviation fitted in TCGA",
        "score_means_tcga": [float(value) for value in score_means],
        "score_scales_tcga": [float(value) for value in score_scales],
        "total_explained_variance_ratio": float(
            model.explained_variance_ratio_.sum()
        ),
        "components": loadings,
    }
    return train_scores, external_scores, metadata
