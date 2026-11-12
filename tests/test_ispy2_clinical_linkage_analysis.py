import pandas as pd
import pytest

from scripts.analyze_ispy2_clinical_linkage import exact_one_to_one_link


def feature_frame(ids=(101, 102)):
    return pd.DataFrame({
        "CLINICAL-TRIAL-SUBJECT-ID": list(ids),
        "VOLUME_TUM_BLU_V10": [10.0] * len(ids),
        "VOLUME_TUM_BLU_V20": [5.0] * len(ids),
        "VOLUME_TUM_BLU_V30": [2.0] * len(ids),
        "VOLUME_TUM_BLU_V40": [1.0] * len(ids),
    })


def clinical_frame(ids=(101, 102)):
    return pd.DataFrame({
        "Patient_ID": list(ids),
        "Arm": ["Paclitaxel"] * len(ids),
        "HR": [0] * len(ids),
        "HER2": [1] * len(ids),
        "MP": [0] * len(ids),
        "pCR": [1] * len(ids),
        "Age_at_Screening": [45] * len(ids),
    })


def test_exact_one_to_one_link_matches_all_feature_rows():
    linked = exact_one_to_one_link(feature_frame(), clinical_frame())
    assert len(linked) == 2
    assert linked["Patient_ID"].tolist() == [101, 102]


def test_link_rejects_duplicate_feature_ids():
    with pytest.raises(ValueError, match="Duplicate participant IDs"):
        exact_one_to_one_link(feature_frame((101, 101)), clinical_frame((101, 102)))


def test_link_rejects_duplicate_clinical_ids():
    with pytest.raises(ValueError, match="Duplicate participant IDs"):
        exact_one_to_one_link(feature_frame((101, 102)), clinical_frame((101, 101)))


def test_link_rejects_unmatched_feature_participant():
    with pytest.raises(ValueError, match="lack clinical matches"):
        exact_one_to_one_link(feature_frame((101, 999)), clinical_frame((101, 102)))
