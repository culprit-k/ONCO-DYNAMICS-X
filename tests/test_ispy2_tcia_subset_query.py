from scripts.query_ispy2_tcia_subset import (
    build_candidate_rows,
    canonical_stage,
    canonical_subject_id,
    extract_tag_values,
    suggested_role,
)


def test_canonical_subject_id_accepts_six_digit_trial_ids_only():
    assert canonical_subject_id(111881) == "111881"
    assert canonical_subject_id("ISPY2_111881") == "111881"
    assert canonical_subject_id("12345") is None
    assert canonical_subject_id("1234567") is None


def test_canonical_stage_preserves_subvisit_separately():
    assert canonical_stage("T0") == "T0"
    assert canonical_stage("T0-TEST") == "T0"
    assert canonical_stage("T1_RETEST") == "T1"
    assert canonical_stage("baseline") is None


def test_extract_tag_values_returns_unique_values():
    payload = [
        {"element": "(0012,0040)", "name": "ClinicalTrialSubjectID", "data": "111881"},
        {"element": "(0012,0050)", "name": "ClinicalTrialTimePointID", "data": "T1"},
        {"element": "(0012,0040)", "name": "ClinicalTrialSubjectID", "data": "111881"},
    ]
    assert extract_tag_values(payload, "(0012,0040)") == {"111881"}
    assert extract_tag_values(payload, "(0012,0050)") == {"T1"}


def test_suggested_roles_do_not_auto_select_series():
    assert suggested_role("ISPY2: VOLSER: uni-lateral cropped: SER", "MR") == "SER_MAP_REVIEW"
    assert suggested_role("ISPY2: VOLSER: uni-lateral cropped: PE1", "MR") == "PE_MAP_REVIEW"
    assert (
        suggested_role("ISPY2: VOLSER: uni-lateral cropped: Analysis Mask", "SEG")
        == "FTV_MASK_REVIEW"
    )
    assert suggested_role("T2 AXIAL", "MR") == "T2_ANATOMICAL_CANDIDATE"
    assert suggested_role("LOCALIZER", "MR") == "MANUAL_SERIES_REVIEW"


def test_candidate_rows_require_exact_dicom_trial_id_and_unambiguous_visit():
    series = [
        {
            "StudyInstanceUID": "study-1",
            "SeriesInstanceUID": "series-1",
            "PatientID": "ISPY2_111881",
            "Modality": "MR",
            "SeriesDescription": "ISPY2: VOLSER: uni-lateral cropped: SER",
            "SeriesNumber": "1000",
            "ImageCount": 10,
        },
        {
            "StudyInstanceUID": "study-2",
            "SeriesInstanceUID": "series-2",
            "PatientID": "ISPY2_222222",
            "Modality": "MR",
            "SeriesDescription": "T2 AXIAL",
        },
        {
            "StudyInstanceUID": "study-3",
            "SeriesInstanceUID": "series-3",
            "PatientID": "ISPY2_111881",
            "Modality": "MR",
            "SeriesDescription": "T2 AXIAL",
        },
    ]
    study_tags = {
        "study-1": {
            "trial_subject_id": "111881",
            "visit_code": "T1",
            "subject_tag_value_count": 1,
            "visit_tag_value_count": 1,
        },
        "study-2": {
            "trial_subject_id": "222222",
            "visit_code": "T0",
            "subject_tag_value_count": 1,
            "visit_tag_value_count": 1,
        },
        "study-3": {
            "trial_subject_id": "111881",
            "visit_code": "T0-TEST",
            "subject_tag_value_count": 1,
            "visit_tag_value_count": 2,
        },
    }
    frame, summary = build_candidate_rows(series, study_tags, {"111881"}, "ISPY2")
    assert frame["series_instance_uid"].tolist() == ["series-1"]
    assert frame.loc[0, "visit_stage"] == "T1"
    assert frame.loc[0, "include"] == 0
    assert summary["matched_subject_count"] == 1
    assert summary["matched_study_count"] == 1
    assert summary["blocked_study_count_missing_or_ambiguous_trial_tags"] == 1
    assert summary["unmatched_target_subject_count"] == 0
