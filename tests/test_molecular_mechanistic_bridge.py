import numpy as np
import pytest


def _bridge():
    from onco_dynamics_x.data_engine.molecular_mechanistic_bridge import (
        MOLECULAR_MODULES,
        describe_mechanistic_link_policy,
        score_molecular_modules,
    )

    return MOLECULAR_MODULES, describe_mechanistic_link_policy, score_molecular_modules


def test_module_scores_are_means_of_prespecified_rank_normalised_features():
    modules, _, score_molecular_modules = _bridge()
    genes = tuple(gene for group in modules.values() for gene in group)
    matrix = np.arange(4 * len(genes), dtype=float).reshape(4, len(genes))
    result = score_molecular_modules(matrix, genes)

    assert result.sample_scores.shape == (4, len(modules))
    assert result.module_names == tuple(modules)
    assert np.isfinite(result.sample_scores).all()
    for column, module_name in enumerate(result.module_names):
        indices = [genes.index(gene) for gene in modules[module_name]]
        assert np.allclose(result.sample_scores[:, column], matrix[:, indices].mean(axis=1))


def test_module_scoring_rejects_insufficient_gene_coverage():
    modules, _, score_molecular_modules = _bridge()
    genes = tuple(modules["proliferation_expression"])
    matrix = np.ones((5, len(genes)))
    with pytest.raises(ValueError, match="gene coverage"):
        score_molecular_modules(matrix, genes)


def test_mechanistic_policy_is_hypothesis_only_and_forbids_direct_parameter_mapping():
    _, describe_mechanistic_link_policy, _ = _bridge()
    policy = describe_mechanistic_link_policy()

    assert policy["status"] == "hypothesis registry only; no parameter updates"
    assert len(policy["module_hypotheses"]) == 4
    assert any("diffusion coefficient" in item for item in policy["prohibited_direct_mappings"])
    assert any("drug sensitivity" in item for item in policy["prohibited_direct_mappings"])
    assert len(policy["required_calibration_sequence"]) == 6


def test_module_scoring_rejects_nonfinite_and_duplicate_inputs():
    modules, _, score_molecular_modules = _bridge()
    genes = tuple(gene for group in modules.values() for gene in group)
    matrix = np.ones((3, len(genes)))
    matrix[0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        score_molecular_modules(matrix, genes)
    with pytest.raises(ValueError, match="unique"):
        score_molecular_modules(np.ones((3, 2)), ("ESR1", "ESR1"))
