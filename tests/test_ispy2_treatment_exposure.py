import pandas as pd
import pytest

from scripts.audit_ispy2_treatment_exposure import align_exposures


def events(time_basis="TCIA_MODIFIED_LONGITUDINAL"):
    return pd.DataFrame(
        {
            "patient_id": ["123456", "123456"],
            "agent": ["Paclitaxel", "Paclitaxel"],
            "administration_datetime": ["1960-01-02T08:00:00", "1960-01-09T08:00:00"],
            "dose": [80.0, 80.0],
            "dose_unit": ["mg/m2", "mg/m2"],
            "route": ["IV", "IV"],
            "administration_status": ["administered", "administered"],
            "time_basis": [time_basis, time_basis],
            "source_reference": ["trial_record_A", "trial_record_B"],
        }
    )


def visits():
    return pd.DataFrame(
        {
            "patient_id": ["123456"],
            "visit_code": ["T2"],
            "scan_datetime": ["1960-01-12T08:00:00"],
            "longitudinal_date_flag": ["MODIFIED"],
        }
    )


def test_aligns_events_only_on_tcia_modified_date_axis():
    aligned, summary = align_exposures(events(), visits())
    assert len(aligned) == 1
    assert aligned.iloc[0]["days_since_last_administration"] == pytest.approx(3.0)
    assert aligned.iloc[0]["cumulative_administered_dose_before_scan"] == pytest.approx(160.0)
    assert summary["time_basis"] == "TCIA_MODIFIED_LONGITUDINAL"


def test_refuses_original_calendar_dates_without_shift_mapping():
    with pytest.raises(ValueError, match="Original-calendar dates cannot be aligned"):
        align_exposures(events("ORIGINAL_CALENDAR"), visits())


def test_requires_modified_date_flag():
    bad = visits()
    bad["longitudinal_date_flag"] = ""
    with pytest.raises(ValueError, match="missing or inconsistent"):
        align_exposures(events(), bad)


def test_excludes_planned_or_cancelled_administrations():
    e = events()
    e.loc[1, "administration_status"] = "planned"
    aligned, summary = align_exposures(e, visits())
    assert aligned.iloc[0]["cumulative_administered_dose_before_scan"] == pytest.approx(80.0)
    assert summary["excluded_non_administered_status_counts"] == {"PLANNED": 1}
