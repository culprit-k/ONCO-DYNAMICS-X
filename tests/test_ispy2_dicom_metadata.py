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
        (0x0117, 0x1012): "FLOAT",
        (0x0117, 0x1014): "SER_min",
        (0x0117, 0x1018): 0,
    })
    ds = FakeDataset({(0x0117, 0x1010): [item]})
    assert extract_series_metadata(ds)["ser_analysis_parameters"]["SER_min"] == 0.0


def test_acquisition_sequence_and_contrast_fields_are_extracted():
    ds = FakeDataset(
        {
            (0x0012, 0x0040): "123456",
            (0x0018, 0x0080): 4.2,
            (0x0018, 0x0081): 1.7,
            (0x0018, 0x0082): 0.8,
            (0x0018, 0x0087): 3.0,
            (0x0018, 0x1314): 12.0,
            (0x0018, 0x0010): "Gadolinium",
            (0x0018, 0x1041): 10.0,
            (0x0018, 0x1042): "101500",
        }
    )
    row = extract_series_metadata(ds)
    assert row["repetition_time_ms"] == 4.2
    assert row["echo_time_ms"] == 1.7
    assert row["inversion_time_ms"] == 0.8
    assert row["magnetic_field_strength_t"] == 3.0
    assert row["flip_angle_degrees"] == 12.0
    assert row["contrast_bolus_agent"] == "Gadolinium"
    assert row["contrast_bolus_volume_ml"] == 10.0
    assert row["contrast_bolus_start_time"] == "101500"


def test_typed_ser_parameters_preserve_float_integer_and_string_values():
    float_item = FakeDataset(
        {
            (0x0117, 0x1012): "FLOAT",
            (0x0117, 0x1014): "PE_threshold",
            (0x0117, 0x1018): "70.5",
        }
    )
    integer_item = FakeDataset(
        {
            (0x0117, 0x1012): "INTEGER",
            (0x0117, 0x1014): "phase_count",
            (0x0117, 0x1019): "12",
        }
    )
    string_item = FakeDataset(
        {
            (0x0117, 0x1012): "STRING",
            (0x0117, 0x1014): "method",
            (0x0117, 0x101A): "AUTO",
        }
    )
    ds = FakeDataset(
        {(0x0117, 0x1010): [float_item, integer_item, string_item]}
    )
    params = extract_series_metadata(ds)["ser_analysis_parameters"]
    assert params["PE_threshold"] == 70.5
    assert params["phase_count"] == 12
    assert params["method"] == "AUTO"


def test_ftv_private_results_preserve_numeric_types():
    ftv_item = FakeDataset(
        {
            (0x0117, 0x10B1): "0.3",
            (0x0117, 0x10B2): "1.8",
            (0x0117, 0x10B3): "125",
            (0x0117, 0x10B4): "1.25",
            (0x0117, 0x10B5): "SER-restricted",
        }
    )
    ds = FakeDataset({(0x0117, 0x10B0): [ftv_item]})
    result = extract_series_metadata(ds)["ftv_results"][0]
    assert result["ser_minimum"] == 0.3
    assert result["ser_maximum"] == 1.8
    assert result["voxel_count"] == 125
    assert result["volume_cc"] == 1.25
    assert result["label"] == "SER-restricted"
