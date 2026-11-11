import numpy as np
import pytest


def _evaluator():
    pytest.importorskip("sklearn")
    from onco_dynamics_x.calibration.pathology import (
        evaluate_module_pathology_association,
    )

    return evaluate_module_pathology_association


def test_pathology_association_uses_patient_level_out_of_fold_predictions():
    evaluate = _evaluator()
    rng = np.random.default_rng(2026)
    labels = ["Negative", "Positive", "Equivocal"]
    features = []
    patient_ids = []
    labels_by_patient = {}
    for index, label in enumerate(labels):
        centre = np.zeros(4)
        centre[index] = 3.0
        for row in range(30):
            patient_id = f"patient-{label}-{row}"
            patient_ids.append(patient_id)
            labels_by_patient[patient_id] = label
            features.append(centre + rng.normal(0.0, 0.3, 4))
    report = evaluate(
        np.asarray(features),
        patient_ids,
        labels_by_patient,
        module_names=("proliferation", "luminal", "her2", "basal"),
        n_splits=5,
        n_bootstrap=30,
        seed=17,
    )

    assert report["patients_evaluated"] == 90
    assert report["metrics"]["accuracy"] > 0.95
    assert report["metrics"]["balanced_accuracy"] > 0.95
    assert report["metrics"]["macro_f1"] > 0.95
    assert report["confusion_matrix"]["class_order"] == sorted(labels)
    assert report["data_governance"]["patient_ids_exported"] is False
    assert "patient-Negative-0" not in str(report)


def test_pathology_association_rejects_duplicate_patient_ids_and_rare_classes():
    evaluate = _evaluator()
    features = np.ones((6, 2))
    labels = {f"P{i}": "A" if i < 3 else "B" for i in range(6)}
    with pytest.raises(ValueError, match="unique"):
        evaluate(
            features,
            ["P0", "P0", "P2", "P3", "P4", "P5"],
            labels,
            module_names=("module_a", "module_b"),
        )
    rare_labels = {f"P{i}": "A" for i in range(5)}
    rare_labels["P5"] = "B"
    with pytest.raises(ValueError, match="at least n_splits"):
        evaluate(
            features,
            [f"P{i}" for i in range(6)],
            rare_labels,
            module_names=("module_a", "module_b"),
            n_splits=3,
        )
