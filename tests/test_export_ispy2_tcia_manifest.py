import pandas as pd
import pytest

from scripts.export_ispy2_tcia_manifest import (
    DEFAULT_REVIEW_NOTE,
    export_reviewed_manifest,
)


def _candidate_rows():
    return pd.DataFrame(
        [
            {
                "trial_subject_id": "111881",
                "visit_code": "T0",
                "visit_stage": "T0",
                "study_instance_uid": "study-0",
                "series_instance_uid": "series-0",
                "modality": "MR",
                "series_description": "T2 AXIAL",
                "include": 1,
                "review_note": "Selected as the axial anatomical reference after image review",
            },
            {
                "trial_subject_id": "111881",
                "visit_code": "T1",
                "visit_stage": "T1",
                "study_instance_uid": "study-1",
                "series_instance_uid": "series-1",
                "modality": "MR",
                "series_description": "VOLSER SER",
                "include": 0,
                "review_note": DEFAULT_REVIEW_NOTE,
            },
        ]
    )


def test_manifest_export_includes_only_human_reviewed_rows():
    manifest, summary = export_reviewed_manifest(_candidate_rows())
    assert manifest.to_dict("list") == {"SeriesInstanceUID": ["series-0"]}
    assert summary["selected_patient_count"] == 1
    assert summary["selected_visit_stage_counts"] == {"T0": 1}


def test_manifest_export_rejects_default_review_note():
    rows = _candidate_rows()
    rows.loc[0, "review_note"] = DEFAULT_REVIEW_NOTE
    with pytest.raises(ValueError, match="human review rationale"):
        export_reviewed_manifest(rows)


def test_manifest_export_rejects_unknown_visit_stage():
    rows = _candidate_rows()
    rows.loc[0, "visit_stage"] = "UNKNOWN"
    with pytest.raises(ValueError, match="verified T0/T1/T2/T3"):
        export_reviewed_manifest(rows)


def test_manifest_export_rejects_duplicate_series_uid():
    rows = _candidate_rows()
    duplicate = rows.iloc[[0]].copy()
    duplicate["series_instance_uid"] = "series-0"
    rows = pd.concat([rows, duplicate], ignore_index=True)
    with pytest.raises(ValueError, match="duplicate SeriesInstanceUID"):
        export_reviewed_manifest(rows)
