import numpy as np
import pytest

from onco_dynamics_x.data_engine.expression import ExpressionMatrix, ExpressionPreprocessor


def test_expression_preprocessor_preserves_gene_order_and_standardizes():
    matrix = ExpressionMatrix(
        sample_ids=("S1", "S2", "S3"),
        gene_ids=("G1", "G2", "G3"),
        values=np.array([[1.0, 5.0, 2.0], [2.0, 5.0, 4.0], [3.0, 5.0, 8.0]]),
    )
    pipeline = ExpressionPreprocessor()
    transformed = pipeline.fit_transform(matrix)
    assert pipeline.selected_gene_ids == ("G1", "G3")
    assert transformed.shape == (3, 2)
    assert np.allclose(transformed.mean(axis=0), 0.0)
    assert np.allclose(transformed.std(axis=0), 1.0)


def test_expression_preprocessor_rejects_missing_training_genes():
    train = ExpressionMatrix(("S1", "S2"), ("G1", "G2"), np.array([[1.0, 2.0], [2.0, 5.0]]))
    test = ExpressionMatrix(("S3",), ("G1",), np.array([[3.0]]))
    pipeline = ExpressionPreprocessor().fit(train)
    with pytest.raises(ValueError, match="missing required genes"):
        pipeline.transform(test)
