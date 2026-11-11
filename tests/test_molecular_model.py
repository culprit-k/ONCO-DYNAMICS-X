import numpy as np
import pytest


def _molecular_model():
    pytest.importorskip("scipy")
    pytest.importorskip("sklearn")
    from onco_dynamics_x.data_engine.molecular_model import (
        fit_tcga_pca_projection,
        rank_normalise_expression,
    )

    return rank_normalise_expression, fit_tcga_pca_projection


def test_rank_normalisation_is_invariant_to_monotonic_gene_scale():
    rank_normalise_expression, _ = _molecular_model()
    generator = np.random.default_rng(12)
    values = generator.lognormal(mean=2.0, sigma=1.2, size=(80, 4))
    transformed = rank_normalise_expression(
        values, ("ESR1", "PGR", "ERBB2", "MKI67")
    )
    rescaled = rank_normalise_expression(
        np.log1p(values), ("ESR1", "PGR", "ERBB2", "MKI67")
    )

    assert transformed.values.shape == values.shape
    assert np.allclose(transformed.values, rescaled.values)
    assert np.all(np.isfinite(transformed.values))


def test_rank_normalisation_uses_median_normal_score_for_sparse_missing_values():
    rank_normalise_expression, _ = _molecular_model()
    values = np.array(
        [
            [1.0, 5.0],
            [2.0, 5.0],
            [3.0, np.nan],
            [4.0, 5.0],
            [5.0, 5.0],
            [6.0, 5.0],
            [7.0, 5.0],
            [8.0, 5.0],
            [9.0, 5.0],
            [10.0, 5.0],
            [11.0, 5.0],
            [12.0, 5.0],
            [13.0, 5.0],
            [14.0, 5.0],
            [15.0, 5.0],
            [16.0, 5.0],
            [17.0, 5.0],
            [18.0, 5.0],
            [19.0, 5.0],
            [20.0, 5.0],
        ]
    )
    result = rank_normalise_expression(
        values, ("ESR1", "PGR"), minimum_observed_fraction=0.95
    )

    assert result.samples_with_imputation == 1
    assert result.missing_fraction_by_gene["ESR1"] == 0.0
    assert result.missing_fraction_by_gene["PGR"] == pytest.approx(0.05)
    assert result.values[2, 1] == 0.0


def test_rank_normalisation_rejects_insufficient_gene_coverage():
    rank_normalise_expression, _ = _molecular_model()
    values = np.full((20, 2), 1.0)
    values[0:2, 1] = np.nan

    with pytest.raises(ValueError, match="observed in only"):
        rank_normalise_expression(
            values, ("ESR1", "PGR"), minimum_observed_fraction=0.95
        )


def test_tcga_fitted_pca_transforms_external_data_with_frozen_components():
    _, fit_tcga_pca_projection = _molecular_model()
    generator = np.random.default_rng(18)
    latent = generator.normal(size=(100, 4))
    mixing = generator.normal(size=(4, 8))
    training = latent @ mixing + generator.normal(scale=0.05, size=(100, 8))
    external = generator.normal(size=(30, 4)) @ mixing + generator.normal(
        scale=0.05, size=(30, 8)
    )
    genes = tuple(f"GENE_{index}" for index in range(8))

    train_scores, external_scores, metadata = fit_tcga_pca_projection(
        training, external, genes, n_components=4, random_state=42
    )

    assert train_scores.shape == (100, 4)
    assert external_scores.shape == (30, 4)
    assert np.isfinite(train_scores).all()
    assert np.isfinite(external_scores).all()
    assert metadata["fit_cohort"] == "TCGA-BRCA"
    assert metadata["external_transform"] == "frozen TCGA PCA components"
    assert metadata["total_explained_variance_ratio"] > 0.95
    assert metadata["score_standardisation"].startswith("zero mean")
    assert np.allclose(train_scores.mean(axis=0), 0.0, atol=1e-10)
    assert np.allclose(train_scores.std(axis=0), 1.0, atol=1e-10)
    assert len(metadata["components"]) == 4
