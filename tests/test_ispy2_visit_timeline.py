import pandas as pd
import pytest

from scripts.build_ispy2_visit_timeline import build_visit_timeline


def inventory_rows(studies=None, flag="MODIFIED"):
    if studies is None:
        studies = [
            ("T0", "1960-01-01T08:00:00", "study0"),
            ("T1", "1960-01-22T08:00:00", "study1"),
            ("T2", "1960-02-19T08:00:00", "study2"),
            ("T3", "1960-04-15T08:00:00", "study3"),
        ]
    return pd.DataFrame(
        [
            {
                "patient_id": "123456",
                "visit_code": visit,
                "study_instance_uid": study,
                "scan_datetime": when,
                "scan_datetime_source": "StudyDate/StudyTime",
                "date_flags": flag,
                "study_datetime_consistent": True,
            }
            for visit, when, study in studies
        ]
    )


def test_timeline_computes_relative_days_from_shifted_dates():
    timeline, summary = build_visit_timeline(inventory_rows())
    days = timeline.sort_values("visit_code")["elapsed_days_from_T0"].astype(float).tolist()
    assert days == [0.0, 21.0, 49.0, 105.0]
    assert summary["visits_with_stage_intervals"] == 4


def test_timeline_refuses_missing_date_modification_flag():
    timeline, summary = build_visit_timeline(inventory_rows(flag=""))
    assert timeline["elapsed_days_from_T0"].isna().all()
    assert summary["interval_status_counts"]["blocked_missing_or_inconsistent_date_flag"] == 4


def test_timeline_refuses_multiple_studies_for_same_visit():
    rows = inventory_rows()
    extra = rows.iloc[[0]].copy()
    extra["study_instance_uid"] = "second_t0"
    extra["scan_datetime"] = "1960-01-02T08:00:00"
    timeline, _ = build_visit_timeline(pd.concat([rows, extra], ignore_index=True))
    assert timeline["elapsed_days_from_T0"].isna().all()
    assert set(timeline["interval_status"]) == {"blocked_multiple_studies_within_visit"}


def test_timeline_refuses_negative_time_from_baseline():
    rows = inventory_rows()
    rows.loc[1, "scan_datetime"] = "1959-12-30T08:00:00"
    timeline, _ = build_visit_timeline(rows)
    t1 = timeline[timeline["visit_code"] == "T1"].iloc[0]
    assert t1["interval_status"] == "blocked_negative_time_from_T0"
    assert pd.isna(t1["elapsed_days_from_T0"])


def test_timeline_normalizes_t1_subvisit_tag():
    rows = inventory_rows()
    rows.loc[1, "visit_code"] = "T1-2"
    rows.loc[1, "visit_stage"] = "T1"
    timeline, _ = build_visit_timeline(rows)
    t1 = timeline[timeline["visit_code"] == "T1"].iloc[0]
    assert t1["visit_code_raw"] == "T1-2"
    assert t1["elapsed_days_from_T0"] == pytest.approx(21.0)


def test_timeline_date_only_precision_uses_calendar_days():
    rows = inventory_rows()
    rows["study_datetime"] = rows["scan_datetime"].str.slice(0, 10) + "T00:00:00"
    rows["study_time_precision"] = "date_only"
    rows["study_datetime_values_json"] = '[ "one-date" ]'
    rows["scan_datetime"] = rows["study_datetime"]
    timeline, _ = build_visit_timeline(rows)
    assert timeline["elapsed_days_from_T0"].astype(float).tolist() == [0.0, 21.0, 49.0, 105.0]
    assert set(timeline["elapsed_time_precision"]) == {"date"}
