"""Expand the selected I-SPY 2 DCE studies to all sibling series in TCIA.

The input is a private cohort-reviewed table containing one original-DCE
SeriesInstanceUID per participant and T0-T3 visit. For each target TCIA
PatientID, query the public TCIA v4 API once, then retain all series whose
StudyInstanceUID belongs to one of the four selected longitudinal visits.

Outputs are private candidate metadata only. Every candidate defaults to
include=0; this script never creates a download manifest or selects a series.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import json
import re
import time
from collections import Counter
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd

TCIA_QUERY = "https://services.cancerimagingarchive.net/services/v4/TCIA/query"
USER_AGENT = "ONCO-DYNAMICS-X/0.1 (I-SPY2 study-series expansion)"
VISITS = {"T0", "T1", "T2", "T3"}


def _get_json(url: str, timeout: int = 120) -> Any:
    request = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    with urlopen(request, timeout=timeout) as response:
        raw = response.read()
    return json.loads(raw.decode("utf-8-sig"))


def canonical_subject_id(value: Any) -> str | None:
    if value is None or pd.isna(value):
        return None
    digits = re.sub(r"\D", "", str(value).strip())
    if 1 <= len(digits) <= 6:
        return digits.zfill(6)
    return None


def classify_series_role(description: Any, modality: Any) -> str:
    """Assign a conservative review label from TCIA description/modality."""
    text = str(description or "").strip()
    lower = text.lower()
    mod = str(modality or "").strip().upper()
    if "analysis mask" in lower:
        return "FTV_ANALYSIS_MASK_CANDIDATE"
    if "volser" in lower and re.search(r"\bpe\s*\d+\b", lower):
        return "PE_MAP_CANDIDATE"
    if "volser" in lower and re.search(r"\bser\b", lower):
        return "SER_MAP_CANDIDATE"
    if "original dce" in lower:
        return "ORIGINAL_DCE_SOURCE"
    if mod == "SEG":
        return "SEGMENTATION_REVIEW"
    if re.search(r"\bt2(?:w| weighted)?\b", lower) or "t2-weighted" in lower:
        return "T2_ANATOMICAL_CANDIDATE"
    if "dynaviews" in lower and "pre" in lower:
        return "DCE_PRECONTRAST_REFERENCE"
    if "dynaviews" in lower and "post" in lower:
        return "DCE_POSTCONTRAST_REFERENCE"
    return "OTHER_SERIES_REVIEW"


def _unwrap_series(payload: Any) -> list[dict[str, Any]]:
    if isinstance(payload, dict):
        for key in ("data", "series", "result"):
            if isinstance(payload.get(key), list):
                payload = payload[key]
                break
    if not isinstance(payload, list):
        raise ValueError("TCIA v4 getSeries response is not a JSON list")
    return [row for row in payload if isinstance(row, dict)]


def _field(row: dict[str, Any], *names: str) -> Any:
    for name in names:
        if row.get(name) not in (None, ""):
            return row[name]
    return None


def fetch_patient_series(patient_id: str, *, attempts: int = 5, sleep_seconds: float = 0.5) -> list[dict[str, Any]]:
    params = urlencode({"Collection": "ISPY2", "PatientID": patient_id, "format": "JSON"})
    url = f"{TCIA_QUERY}/getSeries?{params}"
    for attempt in range(attempts):
        try:
            return _unwrap_series(_get_json(url))
        except HTTPError as exc:
            if exc.code not in {429, 500, 502, 503, 504} or attempt == attempts - 1:
                raise
            retry_after = exc.headers.get("Retry-After") if exc.headers else None
            try:
                delay = max(float(retry_after), 2.0 ** attempt) if retry_after else 2.0 ** attempt
            except ValueError:
                delay = 2.0 ** attempt
            time.sleep(min(delay, 60.0))
        except (URLError, TimeoutError, json.JSONDecodeError):
            if attempt == attempts - 1:
                raise
            time.sleep(min(2.0 ** attempt, 30.0))
    raise RuntimeError("Unreachable retry state")


def build_study_targets(review: pd.DataFrame) -> tuple[dict[str, dict[str, Any]], dict[str, set[str]]]:
    required = {
        "PatientID", "ClinicalTrialSubjectID", "TimePoint",
        "StudyInstanceUID", "SeriesInstanceUID", "SeriesDescription",
    }
    missing = required - set(review.columns)
    if missing:
        raise ValueError(f"Review CSV missing required columns: {sorted(missing)}")

    data = review.copy()
    data["ClinicalTrialSubjectID"] = data["ClinicalTrialSubjectID"].map(canonical_subject_id)
    data["TimePoint"] = data["TimePoint"].astype("string").str.strip().str.upper()
    data["PatientID"] = data["PatientID"].astype("string").str.strip()
    data["StudyInstanceUID"] = data["StudyInstanceUID"].astype("string").str.strip()
    data["SeriesInstanceUID"] = data["SeriesInstanceUID"].astype("string").str.strip()

    if data[list(required)].isna().any().any():
        raise ValueError("Review CSV has missing required values; resolve before study expansion")
    if not data["TimePoint"].isin(VISITS).all():
        raise ValueError("Review TimePoint values must be T0, T1, T2 or T3")
    if data["SeriesDescription"].astype(str).str.contains("original dce", case=False).mean() < 0.99:
        raise ValueError("Input rows must be original-DCE source series (at least 99% by description)")

    if data["SeriesInstanceUID"].duplicated().any():
        raise ValueError("Input review CSV has duplicate original-DCE SeriesInstanceUID values")
    if data.duplicated(["ClinicalTrialSubjectID", "TimePoint"]).any():
        raise ValueError("Input review CSV has multiple original-DCE series for a participant/visit")
    patient_counts = data.groupby("ClinicalTrialSubjectID")["TimePoint"].nunique()
    if not patient_counts.eq(4).all():
        raise ValueError("Input review CSV must contain exactly one original-DCE source for each T0-T3 visit")

    study_map: dict[str, dict[str, Any]] = {}
    patient_studies: dict[str, set[str]] = {}
    for row in data.to_dict("records"):
        study_uid = str(row["StudyInstanceUID"])
        patient_id = str(row["PatientID"])
        subject_id = str(row["ClinicalTrialSubjectID"])
        visit = str(row["TimePoint"])
        if study_uid in study_map:
            old = study_map[study_uid]
            identity = (old["patient_id"], old["subject_id"], old["visit_stage"])
            current = (patient_id, subject_id, visit)
            if identity != current:
                raise ValueError(f"StudyInstanceUID maps to conflicting patient/visit rows: {study_uid}")
        study_map[study_uid] = {
            "patient_id": patient_id,
            "subject_id": subject_id,
            "visit_stage": visit,
            "source_dce_series_uid": str(row["SeriesInstanceUID"]),
        }
        patient_studies.setdefault(patient_id, set()).add(study_uid)
    return study_map, patient_studies


def expand_study_series(
    review: pd.DataFrame,
    *,
    workers: int = 2,
    cache_dir: Path | None = None,
) -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame]:
    if workers < 1 or workers > 4:
        raise ValueError("workers must be between 1 and 4")
    study_map, patient_studies = build_study_targets(review)
    cache = cache_dir
    if cache is not None:
        cache.mkdir(parents=True, exist_ok=True)

    def load_patient(patient_id: str) -> tuple[str, list[dict[str, Any]], str | None]:
        cache_file = cache / f"{re.sub(r'[^A-Za-z0-9_-]', '_', patient_id)}.json" if cache else None
        if cache_file and cache_file.is_file():
            try:
                rows = _unwrap_series(json.loads(cache_file.read_text(encoding="utf-8")))
                return patient_id, rows, None
            except Exception:
                cache_file.unlink(missing_ok=True)
        try:
            rows = fetch_patient_series(patient_id)
            if cache_file:
                cache_file.write_text(json.dumps(rows, ensure_ascii=False), encoding="utf-8")
            return patient_id, rows, None
        except Exception as exc:
            return patient_id, [], f"{type(exc).__name__}: {exc}"

    patients = sorted(patient_studies)
    results: dict[str, list[dict[str, Any]]] = {}
    errors: list[dict[str, str]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(load_patient, patient_id): patient_id for patient_id in patients}
        for future in concurrent.futures.as_completed(futures):
            patient_id, rows, error = future.result()
            results[patient_id] = rows
            if error:
                errors.append({"PatientID": patient_id, "error": error})

    candidates: list[dict[str, Any]] = []
    seen_uids: set[str] = set()
    for patient_id in patients:
        target_studies = patient_studies[patient_id]
        for row in results.get(patient_id, []):
            study_uid = str(_field(row, "StudyInstanceUID", "studyInstanceUID") or "").strip()
            if study_uid not in target_studies:
                continue
            series_uid = str(_field(row, "SeriesInstanceUID", "seriesInstanceUID") or "").strip()
            if not series_uid:
                continue
            if series_uid in seen_uids:
                raise ValueError(f"TCIA API returned duplicate SeriesInstanceUID: {series_uid}")
            seen_uids.add(series_uid)
            target = study_map[study_uid]
            modality = _field(row, "Modality", "modality")
            description = _field(row, "SeriesDescription", "seriesDescription")
            role = classify_series_role(description, modality)
            candidates.append({
                "trial_subject_id": target["subject_id"],
                "patient_id": patient_id,
                "visit_code": target["visit_stage"],
                "visit_stage": target["visit_stage"],
                "study_instance_uid": study_uid,
                "series_instance_uid": series_uid,
                "source_dce_series_uid": target["source_dce_series_uid"],
                "is_source_dce_series": series_uid == target["source_dce_series_uid"],
                "modality": modality,
                "series_number": _field(row, "SeriesNumber", "seriesNumber"),
                "series_description": description,
                "protocol_name": _field(row, "ProtocolName", "protocolName"),
                "image_count": _field(row, "ImageCount", "imageCount"),
                "manufacturer": _field(row, "Manufacturer", "manufacturer"),
                "manufacturer_model_name": _field(row, "ManufacturerModelName", "manufacturerModelName"),
                "suggested_role": role,
                "include": 0,
                "review_note": "Candidate only; verify geometry/content and study association before download",
            })

    frame = pd.DataFrame(candidates)
    if not frame.empty:
        frame = frame.sort_values(
            ["trial_subject_id", "visit_stage", "suggested_role", "series_number"],
            kind="stable",
        ).reset_index(drop=True)
    expected_source_uids = set(review["SeriesInstanceUID"].astype(str).str.strip())
    observed_source_uids = set(frame.loc[frame["is_source_dce_series"], "series_instance_uid"].astype(str)) if not frame.empty else set()
    missing_source = sorted(expected_source_uids - observed_source_uids)
    summary = {
        "target_participant_count": int(review["ClinicalTrialSubjectID"].map(canonical_subject_id).nunique()),
        "target_visit_study_count": int(len(study_map)),
        "patient_query_count": int(len(patients)),
        "patient_queries_successful": int(len(patients) - len(errors)),
        "patient_query_errors": int(len(errors)),
        "series_in_target_studies": int(len(frame)),
        "source_dce_series_expected": int(len(expected_source_uids)),
        "source_dce_series_returned": int(len(observed_source_uids)),
        "source_dce_series_missing": int(len(missing_source)),
        "missing_source_series_uids": missing_source,
        "role_counts": dict(Counter(frame["suggested_role"])) if not frame.empty else {},
        "visit_stage_counts": dict(Counter(frame["visit_stage"])) if not frame.empty else {},
        "selection_policy": "All candidate series default to include=0; human image and metadata review is mandatory.",
        "privacy": "Candidate CSV contains patient IDs and DICOM UIDs; keep private and do not commit",
    }
    error_frame = pd.DataFrame(errors, columns=["PatientID", "error"])
    return frame, summary, error_frame


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review-csv", required=True, type=Path)
    parser.add_argument("--candidate-csv", required=True, type=Path)
    parser.add_argument("--summary-json", required=True, type=Path)
    parser.add_argument("--errors-csv", required=True, type=Path)
    parser.add_argument("--cache-dir", type=Path, help="Private per-patient metadata cache for safe resumption")
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    review = pd.read_csv(args.review_csv, dtype="string")
    frame, summary, errors = expand_study_series(
        review, workers=args.workers, cache_dir=args.cache_dir
    )
    for path in (args.candidate_csv, args.summary_json, args.errors_csv):
        path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.candidate_csv, index=False)
    args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    errors.to_csv(args.errors_csv, index=False)
    print(json.dumps(summary, indent=2))
    print(f"Candidate CSV: {args.candidate_csv}")
    print(f"Query errors CSV: {args.errors_csv}")


if __name__ == "__main__":
    main()
