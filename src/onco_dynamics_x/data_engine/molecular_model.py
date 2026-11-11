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


INTRINSIC_SUBTYPE_CLASSES: tuple[str, ...] = (
    "Basal",
    "Her2",
    "LumA",
    "LumB",
    "Normal",
)


def canonical_intrinsic_subtype(value: str | None) -> str | None:
    """Map study-specific subtype strings to five shared intrinsic classes.

    Extra METABRIC labels such as claudin-low and NC deliberately return None.
    """
    token = str(value or "").strip().upper().replace(" ", "").replace("-", "_")
    if token.startswith("BRCA_"):
        token = token[5:]
    aliases = {
        "BASAL": "Basal",
        "HER2": "Her2",
        "HER": "Her2",
        "LUMA": "LumA",
        "LUMB": "LumB",
        "NORMAL": "Normal",
    }
    return aliases.get(token)


def evaluate_intrinsic_subtype_transfer(
    training_features: NDArray[np.float64],
    training_patient_ids: list[str] | tuple[str, ...],
    training_labels_by_patient: dict[str, str],
    external_features: NDArray[np.float64],
    external_patient_ids: list[str] | tuple[str, ...],
    external_labels_by_patient: dict[str, str],
    *,
    c_value: float = 0.1,
    n_bootstrap: int = 250,
    seed: int = 42,
) -> dict[str, Any]:
    """Train a fixed logistic classifier in TCGA and test shared labels in METABRIC.

    Features must be generated without subtype labels or clinical outcomes.
    Only five shared classes are used. This benchmarks transfer of existing
    study annotations; it is not an official PAM50/Prosigna classifier.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import (
        accuracy_score,
        balanced_accuracy_score,
        confusion_matrix,
        f1_score,
        precision_recall_fscore_support,
    )

    train_x = np.asarray(training_features, dtype=float)
    test_x = np.asarray(external_features, dtype=float)
    train_ids = tuple(str(value) for value in training_patient_ids)
    test_ids = tuple(str(value) for value in external_patient_ids)
    if train_x.ndim != 2 or test_x.ndim != 2:
        raise ValueError("Subtype feature matrices must be two-dimensional.")
    if train_x.shape[0] != len(train_ids) or test_x.shape[0] != len(test_ids):
        raise ValueError("Patient IDs must align one-to-one with feature rows.")
    if train_x.shape[1] != test_x.shape[1] or train_x.shape[1] < 1:
        raise ValueError("Training and external feature dimensions must match.")
    if not np.isfinite(train_x).all() or not np.isfinite(test_x).all():
        raise ValueError("Subtype features must be finite.")
    if len(set(train_ids)) != len(train_ids) or len(set(test_ids)) != len(test_ids):
        raise ValueError("Patient IDs must be unique within each cohort.")
    if not np.isfinite(c_value) or c_value <= 0:
        raise ValueError("c_value must be finite and positive.")
    if n_bootstrap < 0:
        raise ValueError("n_bootstrap cannot be negative.")

    classes = list(INTRINSIC_SUBTYPE_CLASSES)

    def select_labels(
        patient_ids: tuple[str, ...],
        label_map: dict[str, str],
    ) -> tuple[np.ndarray, np.ndarray, dict[str, int], dict[str, int]]:
        indices: list[int] = []
        labels: list[str] = []
        raw_counts: dict[str, int] = {}
        excluded_counts: dict[str, int] = {}
        for index, patient_id in enumerate(patient_ids):
            raw = str(label_map.get(patient_id, "")).strip()
            if not raw:
                continue
            raw_counts[raw] = raw_counts.get(raw, 0) + 1
            canonical = canonical_intrinsic_subtype(raw)
            if canonical is None:
                excluded_counts[raw] = excluded_counts.get(raw, 0) + 1
                continue
            indices.append(index)
            labels.append(canonical)
        return (
            np.asarray(indices, dtype=int),
            np.asarray(labels, dtype=str),
            dict(sorted(raw_counts.items())),
            dict(sorted(excluded_counts.items())),
        )

    train_indices, y_train, train_raw, train_excluded = select_labels(
        train_ids, training_labels_by_patient
    )
    test_indices, y_test, test_raw, test_excluded = select_labels(
        test_ids, external_labels_by_patient
    )
    if len(train_indices) < 100 or len(test_indices) < 100:
        raise ValueError("Insufficient labelled patients for five-class subtype transfer.")
    train_counts = {
        label: int(np.count_nonzero(y_train == label)) for label in classes
    }
    test_counts = {
        label: int(np.count_nonzero(y_test == label)) for label in classes
    }
    if any(train_counts[label] < 10 for label in classes):
        raise ValueError("At least ten TCGA training patients per shared subtype are required.")
    if any(test_counts[label] < 10 for label in classes):
        raise ValueError("At least ten METABRIC test patients per shared subtype are required.")

    classifier = LogisticRegression(
        C=float(c_value),
        penalty="l2",
        solver="lbfgs",
        class_weight="balanced",
        max_iter=3000,
        random_state=seed,
    )
    classifier.fit(train_x[train_indices], y_train)
    predictions = classifier.predict(test_x[test_indices])
    accuracy = float(accuracy_score(y_test, predictions))
    balanced_accuracy = float(balanced_accuracy_score(y_test, predictions))
    macro_f1 = float(
        f1_score(y_test, predictions, labels=classes, average="macro", zero_division=0)
    )
    precision, recall, class_f1, support = precision_recall_fscore_support(
        y_test, predictions, labels=classes, zero_division=0
    )
    matrix = confusion_matrix(y_test, predictions, labels=classes)

    majority_class = max(train_counts, key=train_counts.get)
    majority_prediction = np.full(len(y_test), majority_class, dtype=f"<U{max(len(label) for label in classes)}")
    majority_accuracy = float(accuracy_score(y_test, majority_prediction))
    majority_macro_f1 = float(
        f1_score(
            y_test, majority_prediction, labels=classes,
            average="macro", zero_division=0,
        )
    )

    rng = np.random.default_rng(seed)
    accuracy_boot: list[float] = []
    balanced_boot: list[float] = []
    macro_f1_boot: list[float] = []
    n_test = len(y_test)
    for _ in range(n_bootstrap):
        sample = rng.integers(0, n_test, size=n_test)
        accuracy_boot.append(float(accuracy_score(y_test[sample], predictions[sample])))
        balanced_boot.append(
            float(balanced_accuracy_score(y_test[sample], predictions[sample]))
        )
        macro_f1_boot.append(
            float(
                f1_score(
                    y_test[sample], predictions[sample], labels=classes,
                    average="macro", zero_division=0,
                )
            )
        )

    def interval(values: list[float]) -> list[float] | None:
        if len(values) < max(30, int(0.5 * n_bootstrap)):
            return None
        return [
            float(np.quantile(values, 0.025)),
            float(np.quantile(values, 0.975)),
        ]

    per_class = {}
    for index, label in enumerate(classes):
        per_class[label] = {
            "precision": float(precision[index]),
            "recall": float(recall[index]),
            "f1": float(class_f1[index]),
            "support": int(support[index]),
        }

    return {
        "task": "five-class intrinsic-subtype label transfer",
        "training_cohort": "TCGA-BRCA",
        "external_evaluation_cohort": "METABRIC",
        "class_order": classes,
        "features": "50-gene panel, cohort-wise outcome-blind rank-normalised values",
        "classifier": {
            "family": "L2 logistic regression",
            "C": float(c_value),
            "class_weight": "balanced",
            "solver": "lbfgs",
            "max_iter": 3000,
            "hyperparameter_tuning": False,
        },
        "training": {
            "patients_with_shared_subtype_labels": int(len(train_indices)),
            "shared_class_counts": train_counts,
            "raw_annotation_value_counts": train_raw,
            "extra_or_nonshared_labels_excluded": train_excluded,
        },
        "external_evaluation": {
            "patients_with_shared_subtype_labels": int(len(test_indices)),
            "shared_class_counts": test_counts,
            "raw_annotation_value_counts": test_raw,
            "extra_or_nonshared_labels_excluded": test_excluded,
            "accuracy": accuracy,
            "accuracy_bootstrap_95ci": interval(accuracy_boot),
            "balanced_accuracy": balanced_accuracy,
            "balanced_accuracy_bootstrap_95ci": interval(balanced_boot),
            "macro_f1": macro_f1,
            "macro_f1_bootstrap_95ci": interval(macro_f1_boot),
            "majority_class_baseline": majority_class,
            "majority_class_baseline_accuracy": majority_accuracy,
            "majority_class_baseline_macro_f1": majority_macro_f1,
            "per_class_metrics": per_class,
            "confusion_matrix": {
                "class_order": classes,
                "counts": matrix.astype(int).tolist(),
            },
        },
        "limitations": [
            "Labels are existing study-provided annotations, not independent biological ground truth.",
            "Features overlap the historical subtype panel, so this measures annotation transfer rather than independent subtype discovery.",
            "Cohort-wise rank normalisation is transductive and not suitable for single-sample deployment.",
            "METABRIC claudin-low and NC labels are excluded from this five-class comparison and counted separately.",
            "This is not an official PAM50/Prosigna classifier and does not establish clinical utility or treatment response.",
        ],
    }
