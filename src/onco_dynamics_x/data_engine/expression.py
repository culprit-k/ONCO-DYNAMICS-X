"""Gene-expression matrix validation and train-only preprocessing."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class ExpressionMatrix:
    """A samples-by-genes expression matrix with explicit sample and gene IDs."""

    sample_ids: tuple[str, ...]
    gene_ids: tuple[str, ...]
    values: FloatArray

    def validate(self) -> None:
        if self.values.ndim != 2:
            raise ValueError("Expression values must be a two-dimensional matrix.")
        if self.values.shape != (len(self.sample_ids), len(self.gene_ids)):
            raise ValueError("Expression matrix shape does not match sample/gene identifiers.")
        if len(set(self.sample_ids)) != len(self.sample_ids):
            raise ValueError("sample_ids must be unique.")
        if len(set(self.gene_ids)) != len(self.gene_ids):
            raise ValueError("gene_ids must be unique.")
        if any(not item for item in self.sample_ids + self.gene_ids):
            raise ValueError("Sample and gene identifiers cannot be empty.")
        if not np.isfinite(self.values).all():
            raise ValueError("Expression values must be finite.")


@dataclass
class ExpressionPreprocessor:
    """Variance filter and z-score scaler fitted only on a training matrix."""

    minimum_variance: float = 1e-8
    selected_gene_ids: tuple[str, ...] = ()
    means: FloatArray | None = None
    scales: FloatArray | None = None

    def fit(self, matrix: ExpressionMatrix) -> "ExpressionPreprocessor":
        matrix.validate()
        if self.minimum_variance < 0:
            raise ValueError("minimum_variance must be non-negative.")
        variances = np.var(matrix.values, axis=0, ddof=0)
        selected = variances > self.minimum_variance
        if not np.any(selected):
            raise ValueError("No genes passed the variance filter.")
        filtered = matrix.values[:, selected]
        self.selected_gene_ids = tuple(
            gene for gene, keep in zip(matrix.gene_ids, selected, strict=True) if keep
        )
        self.means = np.mean(filtered, axis=0)
        raw_scales = np.std(filtered, axis=0, ddof=0)
        self.scales = np.where(raw_scales > 0, raw_scales, 1.0)
        return self

    def transform(self, matrix: ExpressionMatrix) -> FloatArray:
        matrix.validate()
        if self.means is None or self.scales is None or not self.selected_gene_ids:
            raise RuntimeError("ExpressionPreprocessor must be fitted before transform.")
        index_by_gene = {gene: index for index, gene in enumerate(matrix.gene_ids)}
        missing = [gene for gene in self.selected_gene_ids if gene not in index_by_gene]
        if missing:
            raise ValueError(f"Input matrix is missing required genes: {', '.join(missing[:5])}.")
        indices = [index_by_gene[gene] for gene in self.selected_gene_ids]
        filtered = matrix.values[:, indices]
        return (filtered - self.means) / self.scales

    def fit_transform(self, matrix: ExpressionMatrix) -> FloatArray:
        return self.fit(matrix).transform(matrix)
