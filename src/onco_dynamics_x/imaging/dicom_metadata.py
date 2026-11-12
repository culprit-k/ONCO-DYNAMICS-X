"""I-SPY 2 DICOM metadata parsing with explicit date and provenance handling.

This module is intentionally importable without pydicom. It accepts pydicom
Dataset objects or small tag-keyed test doubles.
"""
from __future__ import annotations

import re
from datetime import date, datetime, time, timedelta, timezone
from typing import Any

DICOM_TAGS = {
    "trial_subject_id": (0x0012, 0x0040),
    "visit_code": (0x0012, 0x0050),
    "visit_description": (0x0012, 0x0051),
    "acquisition_datetime": (0x0008, 0x002A),
    "acquisition_date": (0x0008, 0x0022),
    "acquisition_time": (0x0008, 0x0032),
    "series_date": (0x0008, 0x0021),
    "series_time": (0x0008, 0x0031),
    "study_date": (0x0008, 0x0020),
    "study_time": (0x0008, 0x0030),
    "timezone_offset": (0x0008, 0x0201),
    "longitudinal_date_flag": (0x0028, 0x0303),
    "dce_phase_count": (0x0117, 0x1030),
    "dce_phase_duration_s": (0x0117, 0x1031),
    "dce_acquisition_start_times_s": (0x0117, 0x1032),
    "dce_injection_time": (0x0117, 0x1033),
    "dce_effective_delays_s": (0x0117, 0x1034),
    "dce_ser_timing_indices": (0x0117, 0x1035),
    "dce_timing_method": (0x0117, 0x103A),
    "dce_timing_comments": (0x0117, 0x103B),
    "ser_parameter_sequence": (0x0117, 0x1010),
    "ftv_result_sequence": (0x0117, 0x10B0),
}

def _unwrap(value: Any) -> Any:
    return getattr(value, "value", value)


def tag_value(dataset: Any, tag: tuple[int, int], default: Any = None) -> Any:
    """Read a DICOM tag from a pydicom Dataset or a mapping test double."""
    try:
        value = dataset.get(tag, default)
    except (AttributeError, TypeError):
        value = default
    if value is default:
        return default
    return _unwrap(value)


def _text(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        return "\\".join(str(x).strip() for x in value)
    out = str(value).strip()
    return out or None


def _timezone_from_offset(value: Any) -> timezone | None:
    raw = _text(value)
    if not raw:
        return None
    match = re.fullmatch(r"([+-])(\d{2})(\d{2})", raw)
    if not match:
        return None
    minutes = int(match.group(2)) * 60 + int(match.group(3))
    if match.group(1) == "-":
        minutes = -minutes
    return timezone(timedelta(minutes=minutes))


def _parse_da(value: Any) -> date | None:
    raw = _text(value)
    if not raw:
        return None
    digits = re.sub(r"\D", "", raw)
    if len(digits) != 8:
        return None
    try:
        return date(int(digits[:4]), int(digits[4:6]), int(digits[6:8]))
    except ValueError:
        return None


def _parse_tm(value: Any) -> time | None:
    raw = _text(value)
    if not raw:
        return None
    # DICOM TM permits HH, HHMM, HHMMSS and fractional seconds.
    body = raw.split("+")[0].split("-")[0]
    whole, dot, fraction = body.partition(".")
    whole = re.sub(r"\D", "", whole)
    if len(whole) not in (2, 4, 6):
        return None
    whole = whole.ljust(6, "0")
    micros = int((re.sub(r"\D", "", fraction) + "000000")[:6]) if dot else 0
    try:
        return time(int(whole[:2]), int(whole[2:4]), int(whole[4:6]), micros)
    except ValueError:
        return None


def parse_dicom_dt(value: Any, timezone_offset: Any = None) -> datetime | None:
    """Parse DICOM DT with optional UTC offset and fractional seconds."""
    raw = _text(value)
    if not raw:
        return None
    match = re.fullmatch(
        r"(\d{4})(\d{2})(\d{2})(\d{2})?(\d{2})?(\d{2})?(?:\.(\d{1,6}))?([+-]\d{4})?",
        raw,
    )
    if not match:
        return None
    year, month, day, hh, mm, ss, frac, explicit_offset = match.groups()
    zone = _timezone_from_offset(explicit_offset or timezone_offset)
    microsecond = int(((frac or "") + "000000")[:6])
    try:
        return datetime(
            int(year), int(month), int(day),
            int(hh or 0), int(mm or 0), int(ss or 0), microsecond,
            tzinfo=zone,
        )
    except ValueError:
        return None


def combine_dicom_date_time(
    date_value: Any, time_value: Any, timezone_offset: Any = None
) -> tuple[datetime | None, str | None]:
    """Combine DICOM DA/TM; explicitly report if only a date was available."""
    parsed_date = _parse_da(date_value)
    if parsed_date is None:
        return None, None
    parsed_time = _parse_tm(time_value)
    zone = _timezone_from_offset(timezone_offset)
    if parsed_time is None:
        return datetime.combine(parsed_date, time.min, tzinfo=zone), "date_only"
    return datetime.combine(parsed_date, parsed_time, tzinfo=zone), "date_and_time"


def _choose_datetime(dataset: Any) -> tuple[datetime | None, str | None, str | None]:
    offset = tag_value(dataset, DICOM_TAGS["timezone_offset"])
    direct = tag_value(dataset, DICOM_TAGS["acquisition_datetime"])
    parsed = parse_dicom_dt(direct, offset)
    if parsed is not None:
        return parsed, "AcquisitionDateTime (0008,002A)", "date_and_time"

    candidates = [
        ("acquisition_date", "acquisition_time", "AcquisitionDate/AcquisitionTime"),
        ("series_date", "series_time", "SeriesDate/SeriesTime"),
        ("study_date", "study_time", "StudyDate/StudyTime"),
    ]
    for date_key, time_key, label in candidates:
        parsed, precision = combine_dicom_date_time(
            tag_value(dataset, DICOM_TAGS[date_key]),
            tag_value(dataset, DICOM_TAGS[time_key]),
            offset,
        )
        if parsed is not None:
            return parsed, label, precision
    return None, None, None


def extract_subject_id(dataset: Any) -> tuple[str | None, str | None]:
    """Return the six-digit trial ID and the DICOM tag used to obtain it."""
    raw = _text(tag_value(dataset, DICOM_TAGS["trial_subject_id"]))
    if raw:
        digits = re.sub(r"\D", "", raw)
        if 1 <= len(digits) <= 6:
            return digits.zfill(6), "ClinicalTrialSubjectID (0012,0040)"
        trailing_id = re.search(r"(?<!\d)(\d{6})\s*$", raw)
        if trailing_id:
            return trailing_id.group(1), "ClinicalTrialSubjectID (0012,0040)"
        # Do not guess which six-digit block is the trial ID in a longer numeric string.
    patient_id = _text(tag_value(dataset, (0x0010, 0x0020)))
    if patient_id:
        digits = re.sub(r"\D", "", patient_id)
        if 1 <= len(digits) <= 6:
            return digits.zfill(6), "PatientID (0010,0020) fallback"
        trailing_id = re.search(r"(?<!\d)(\d{6})\s*$", patient_id)
        if trailing_id:
            return trailing_id.group(1), "PatientID (0010,0020) fallback"
    return None, None


def extract_series_metadata(dataset: Any) -> dict[str, Any]:
    """Extract the fields needed for a private local series inventory."""
    subject_id, subject_source = extract_subject_id(dataset)
    acquisition_dt, acquisition_source, acquisition_precision = _choose_datetime(dataset)
    study_dt, study_precision = combine_dicom_date_time(
        tag_value(dataset, DICOM_TAGS["study_date"]),
        tag_value(dataset, DICOM_TAGS["study_time"]),
        tag_value(dataset, DICOM_TAGS["timezone_offset"]),
    )
    visit = _text(tag_value(dataset, DICOM_TAGS["visit_code"]))
    if visit:
        visit = visit.upper()
    stage_match = re.match(r"^(T[0-3])(?:[-_].*)?$", visit or "")
    visit_stage = stage_match.group(1) if stage_match else None
    return {
        "patient_id": subject_id,
        "patient_id_source": subject_source,
        "visit_code": visit,
        "visit_stage": visit_stage,
        "visit_description": _text(tag_value(dataset, DICOM_TAGS["visit_description"])),
        "study_instance_uid": _text(tag_value(dataset, (0x0020, 0x000D))),
        "series_instance_uid": _text(tag_value(dataset, (0x0020, 0x000E))),
        "sop_instance_uid": _text(tag_value(dataset, (0x0008, 0x0018))),
        "modality": _text(tag_value(dataset, (0x0008, 0x0060))),
        "series_description": _text(tag_value(dataset, (0x0008, 0x103E))),
        "protocol_name": _text(tag_value(dataset, (0x0018, 0x1030))),
        "series_number": _text(tag_value(dataset, (0x0020, 0x0011))),
        "study_datetime": study_dt.isoformat() if study_dt else None,
        "study_time_precision": study_precision,
        "acquisition_datetime": acquisition_dt.isoformat() if acquisition_dt else None,
        "acquisition_datetime_source": acquisition_source,
        "acquisition_time_precision": acquisition_precision,
        "longitudinal_date_flag": _text(
            tag_value(dataset, DICOM_TAGS["longitudinal_date_flag"])
        ),
        "frame_of_reference_uid": _text(tag_value(dataset, (0x0020, 0x0052))),
        "rows": tag_value(dataset, (0x0028, 0x0010)),
        "columns": tag_value(dataset, (0x0028, 0x0011)),
        "number_of_frames": tag_value(dataset, (0x0028, 0x0008)),
        "pixel_spacing": _text(tag_value(dataset, (0x0028, 0x0030))),
        "rescale_intercept": tag_value(dataset, (0x0028, 0x1052)),
        "rescale_slope": tag_value(dataset, (0x0028, 0x1053)),
        "slice_thickness_mm": tag_value(dataset, (0x0018, 0x0050)),
        "spacing_between_slices_mm": tag_value(dataset, (0x0018, 0x0088)),
        "dce_phase_count": tag_value(dataset, DICOM_TAGS["dce_phase_count"]),
        "dce_phase_duration_s": tag_value(dataset, DICOM_TAGS["dce_phase_duration_s"]),
        "dce_acquisition_start_times_s": _text(
            tag_value(dataset, DICOM_TAGS["dce_acquisition_start_times_s"])
        ),
        "dce_injection_time": _text(tag_value(dataset, DICOM_TAGS["dce_injection_time"])),
        "dce_effective_delays_s": _text(
            tag_value(dataset, DICOM_TAGS["dce_effective_delays_s"])
        ),
        "dce_ser_timing_indices": _text(
            tag_value(dataset, DICOM_TAGS["dce_ser_timing_indices"])
        ),
        "dce_timing_method": _text(tag_value(dataset, DICOM_TAGS["dce_timing_method"])),
        "dce_timing_comments": _text(tag_value(dataset, DICOM_TAGS["dce_timing_comments"])),
        "ftv_results": extract_ftv_results(dataset),
        "ser_analysis_parameters": extract_ser_parameters(dataset),
    }


def _sequence(dataset: Any, tag: tuple[int, int]) -> list[Any]:
    value = tag_value(dataset, tag, [])
    if value is None:
        return []
    return list(value) if isinstance(value, (list, tuple)) else []


def extract_ftv_results(dataset: Any) -> list[dict[str, Any]]:
    """Read the documented UCSF FTV sequence and retain both FTV thresholds."""
    tag = DICOM_TAGS["ftv_result_sequence"]
    rows = []
    for item in _sequence(dataset, tag):
        rows.append({
            "ser_minimum": _text(tag_value(item, (0x0117, 0x10B1))),
            "ser_maximum": _text(tag_value(item, (0x0117, 0x10B2))),
            "voxel_count": _text(tag_value(item, (0x0117, 0x10B3))),
            "volume_cc": _text(tag_value(item, (0x0117, 0x10B4))),
            "label": _text(tag_value(item, (0x0117, 0x10B5))),
        })
    return rows


def extract_ser_parameters(dataset: Any) -> dict[str, str]:
    """Parse named VOLSER parameters without assuming a default PE threshold."""
    tag = DICOM_TAGS["ser_parameter_sequence"]
    values: dict[str, str] = {}
    for item in _sequence(dataset, tag):
        name = _text(tag_value(item, (0x0117, 0x1014)))
        if not name:
            continue
        value = None
        for value_tag in ((0x0117, 0x1018), (0x0117, 0x1019), (0x0117, 0x101A)):
            candidate = tag_value(item, value_tag)
            if candidate is not None:
                value = candidate
                break
        values[name] = _text(value) or ""
    return values


def elapsed_days(start: str, end: str) -> float:
    """Difference between parsed ISO datetimes; expects same patient date offset."""
    left = datetime.fromisoformat(start)
    right = datetime.fromisoformat(end)
    if (left.tzinfo is None) != (right.tzinfo is None):
        raise ValueError("Cannot compare naive and timezone-aware timestamps")
    if left.tzinfo and right.tzinfo and left.utcoffset() != right.utcoffset():
        right = right.astimezone(left.tzinfo)
    return (right - left).total_seconds() / 86400.0
