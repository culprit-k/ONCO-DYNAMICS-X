from datetime import datetime

from onco_dynamics_x.imaging.dicom_metadata import (
    combine_dicom_date_time,
    elapsed_days,
    extract_series_metadata,
    extract_subject_id,
    parse_dicom_dt,
)


class FakeDataset(dict):
    pass


def test_trial_patient_id_and_visit_code_are_read_from_documented_tags():
    ds = FakeDataset({
        (0x0012, 0x0040): "123456",
        (0x0012, 0x0050): "T2",
        (0x0012, 0x0051): "inter-regimen",
    })
    assert extract_subject_id(ds) == ("123456", "ClinicalTrialSubjectID (0012,0040)")
    row = extract_series_metadata(ds)
    assert row["visit_code"] == "T2"
    assert row["visit_description"] == "inter-regimen"


def test_patient_id_falls_back_to_patient_id_tag():
    ds = FakeDataset({(0x0010, 0x0020): "ACRIN-6698-004321"})
    assert extract_subject_id(ds)[0] == "004321"


def test_dicom_dt_preserves_offset_and_fractional_seconds():
    parsed = parse_dicom_dt("20251025103045.125+0500")
    assert parsed.isoformat() == "2025-10-25T10:30:45.125000+05:00"


def test_dicom_da_tm_and_date_only_precision():
    parsed, precision = combine_dicom_date_time("20251025", "103045", "+0500")
    assert parsed.isoformat() == "2025-10-25T10:30:45+05:00"
    assert precision == "date_and_time"
    date_only, precision = combine_dicom_date_time("20251025", None, None)
    assert date_only == datetime(2025, 10, 25)
    assert precision == "date_only"


def test_elapsed_days_between_shifted_longitudinal_dates():
    assert elapsed_days("1960-01-02T12:00:00", "1960-01-12T12:00:00") == 10.0


def test_visit_code_suffix_is_preserved_and_normalized_to_stage():
    ds = FakeDataset({(0x0012, 0x0040): "123456", (0x0012, 0x0050): "T1-2"})
    row = extract_series_metadata(ds)
    assert row["visit_code"] == "T1-2"
    assert row["visit_stage"] == "T1"


def test_ambiguous_long_trial_identifier_is_not_guessed():
    ds = FakeDataset({(0x0012, 0x0040): "1234567890"})
    assert extract_subject_id(ds) == (None, None)


def test_zero_value_ser_parameter_is_preserved():
    item = FakeDataset({
        (0x0117, 0x1014): "SER_min",
        (0x0117, 0x1018): 0,
    })
    ds = FakeDataset({(0x0117, 0x1010): [item]})
    assert extract_series_metadata(ds)["ser_analysis_parameters"]["SER_min"] == "0"
