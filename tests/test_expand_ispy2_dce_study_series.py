import pandas as pd
import pytest

from scripts.expand_ispy2_dce_study_series import (
    build_study_targets,
    classify_series_role,
    expand_study_series,
)


def _review():
    rows = []
    for visit, study, series in [
        ("T0", "study-0", "source-0"),
        ("T1", "study-1", "source-1"),
        ("T2", "study-2", "source-2"),
        ("T3", "study-3", "source-3"),
    ]:
        rows.append({
            "PatientID": "ISPY2-111881",
            "ClinicalTrialSubjectID": "111881",
            "TimePoint": visit,
            "StudyInstanceUID": study,
            "SeriesInstanceUID": series,
            "SeriesDescription": "ISPY2: VOLSER: uni-lateral cropped: original DCE",
        })
    return pd.DataFrame(rows)


def test_role_classifier_identifies_derived_dce_objects():
    assert classify_series_role("ISPY2: VOLSER: cropped: PE1", "MR") == "PE_MAP_CANDIDATE"
    assert classify_series_role("ISPY2: VOLSER: cropped: SER", "MR") == "SER_MAP_CANDIDATE"
    assert classify_series_role("ISPY2: Analysis Mask", "MR") == "FTV_ANALYSIS_MASK_CANDIDATE"
    assert classify_series_role("ISPY2: VOLSER: original DCE", "MR") == "ORIGINAL_DCE_SOURCE"
    assert classify_series_role("T2 Axial", "MR") == "T2_ANATOMICAL_CANDIDATE"


def test_study_targets_require_one_original_dce_per_visit():
    targets, by_patient = build_study_targets(_review())
    assert len(targets) == 4
    assert by_patient == {"ISPY2-111881": {"study-0", "study-1", "study-2", "study-3"}}


def test_study_targets_reject_duplicate_visit_rows():
    review = pd.concat([_review(), _review().iloc[[0]].assign(
        StudyInstanceUID="study-0-retest", SeriesInstanceUID="source-0-retest"
    )], ignore_index=True)
    with pytest.raises(ValueError, match="multiple original-DCE series"):
        build_study_targets(review)


def test_expand_keeps_only_series_from_target_studies(monkeypatch):
    responses = {
        "ISPY2-111881": [
            {
                "StudyInstanceUID": "study-0",
                "SeriesInstanceUID": "source-0",
                "SeriesDescription": "ISPY2: VOLSER: original DCE",
                "Modality": "MR",
                "ImageCount": 1232,
            },
            {
                "StudyInstanceUID": "study-0",
                "SeriesInstanceUID": "ser-0",
                "SeriesDescription": "ISPY2: VOLSER: SER",
                "Modality": "MR",
                "ImageCount": 176,
            },
            {
                "StudyInstanceUID": "unselected-study",
                "SeriesInstanceUID": "unselected-series",
                "SeriesDescription": "T2 AXIAL",
                "Modality": "MR",
            },
            *[
                {
                    "StudyInstanceUID": f"study-{i}",
                    "SeriesInstanceUID": f"source-{i}",
                    "SeriesDescription": "ISPY2: VOLSER: original DCE",
                    "Modality": "MR",
                }
                for i in (1, 2, 3)
            ],
        ]
    }

    def fake_fetch(patient_id, **kwargs):
        return responses[patient_id]

    monkeypatch.setattr(
        "scripts.expand_ispy2_dce_study_series.fetch_patient_series", fake_fetch
    )
    candidates, summary, errors = expand_study_series(_review(), workers=1)

    assert len(candidates) == 5
    assert set(candidates["series_instance_uid"]) == {
        "source-0", "ser-0", "source-1", "source-2", "source-3"
    }
    assert summary["source_dce_series_missing"] == 0
    assert summary["target_participant_count"] == 1
    assert summary["target_visit_study_count"] == 4
    assert set(candidates["include"]) == {0}
    assert errors.empty
