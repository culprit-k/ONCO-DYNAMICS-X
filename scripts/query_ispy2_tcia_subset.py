"""Query a patient-selected I-SPY 2 series shortlist from public TCIA metadata.

The script uses the public TCIA v4 series API, then verifies the clinical-trial
subject and visit tags for one representative series in each StudyInstanceUID
through the NBIA guest-token API. It writes patient-level candidates locally;
never commit its CSV output or token.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import json
import re
import time
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd

TCIA_V4 = "https://services.cancerimagingarchive.net/services/v4/TCIA/query"
NBIA_BASE = "https://services.cancerimagingarchive.net/nbia-api"
USER_AGENT = "ONCO-DYNAMICS-X/0.1 (research metadata query)"
SUBJECT_TAG = "(0012,0040)"
VISIT_TAG = "(0012,0050)"
STAGE_PATTERN = re.compile(r"^(T[0-3])(?:[-_].*)?$", re.IGNORECASE)


def _request_bytes(url: str, *, data: dict[str, str] | None = None,
                   headers: dict[str, str] | None = None, timeout: int = 60) -> bytes:
    request_headers = {"User-Agent": USER_AGENT}
    if headers:
        request_headers.update(headers)
    body = urlencode(data).encode("utf-8") if data is not None else None
    request = Request(url, data=body, headers=request_headers)
    with urlopen(request, timeout=timeout) as response:
        return response.read()


def _request_json(url: str, *, data: dict[str, str] | None = None,
                  headers: dict[str, str] | None = None, timeout: int = 60) -> Any:
    raw = _request_bytes(url, data=data, headers=headers, timeout=timeout)
    try:
        return json.loads(raw.decode("utf-8-sig"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        preview = raw[:300].decode("utf-8", errors="replace")
        raise ValueError(f"TCIA returned a non-JSON response for {url}: {preview}") from exc


def _request_with_retry(url: str, *, headers: dict[str, str], workersafe: bool = True,
                        attempts: int = 5) -> Any:
    for attempt in range(attempts):
        try:
            return _request_json(url, headers=headers)
        except HTTPError as exc:
            if exc.code not in {429, 500, 502, 503, 504} or attempt == attempts - 1:
                raise
            time.sleep(min(60.0, 2.0 ** attempt))
        except (URLError, TimeoutError):
            if attempt == attempts - 1:
                raise
            time.sleep(min(30.0, 2.0 ** attempt))
    raise RuntimeError("Unreachable retry state")


def fetch_public_series(collection: str) -> list[dict[str, Any]]:
    url = f"{TCIA_V4}/getSeries?{urlencode({'Collection': collection, 'format': 'JSON'})}"
    payload = _request_json(url, timeout=180)
    if isinstance(payload, dict):
        for key in ("data", "series", "result"):
            if isinstance(payload.get(key), list):
                payload = payload[key]
                break
    if not isinstance(payload, list):
        raise ValueError("Unexpected getSeries response: expected a JSON list of series")
    return [row for row in payload if isinstance(row, dict)]


def request_guest_token() -> str:
    payload = _request_json(
        f"{NBIA_BASE}/oauth/token",
        data={
            "username": "nbia_guest",
            "password": "",
            "client_id": "NBIA",
            "grant_type": "password",
        },
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        timeout=60,
    )
    token = payload.get("access_token") if isinstance(payload, dict) else None
    if not token:
        raise RuntimeError("TCIA did not return a public NBIA guest access token")
    return str(token)


def _normalise_tag(tag: Any) -> str:
    return re.sub(r"\s+", "", str(tag or "")).upper()


def extract_tag_values(payload: Any, wanted_tag: str) -> set[str]:
    """Return distinct non-empty values for a DICOM tag from getDicomTags JSON."""
    if isinstance(payload, dict):
        rows = payload.get("data", payload.get("tags", [payload]))
    else:
        rows = payload
    if not isinstance(rows, list):
        return set()
    wanted = _normalise_tag(wanted_tag)
    values: set[str] = set()
    for row in rows:
        if not isinstance(row, dict):
            continue
        element = row.get("element", row.get("tag", row.get("Tag", "")))
        if _normalise_tag(element) == wanted:
            value = row.get("data", row.get("value", row.get("Value")))
            if value is not None and str(value).strip():
                values.add(str(value).strip())
    return values


def fetch_study_trial_tags(series_uid: str, token: str) -> dict[str, Any]:
    url = f"{NBIA_BASE}/services/getDicomTags?{urlencode({'SeriesUID': series_uid})}"
    payload = _request_with_retry(url, headers={"Authorization": f"Bearer {token}"})
    subject_values = extract_tag_values(payload, SUBJECT_TAG)
    visit_values = extract_tag_values(payload, VISIT_TAG)
    return {
        "trial_subject_id": next(iter(subject_values)) if len(subject_values) == 1 else None,
        "visit_code": next(iter(visit_values)) if len(visit_values) == 1 else None,
        "subject_tag_value_count": len(subject_values),
        "visit_tag_value_count": len(visit_values),
    }


def canonical_subject_id(value: Any) -> str | None:
    if value is None:
        return None
    match = re.search(r"(?<!\d)(\d{6})$", str(value).strip())
    return match.group(1) if match else None


def canonical_stage(visit_code: Any) -> str | None:
    match = STAGE_PATTERN.match(str(visit_code or "").strip())
    return match.group(1).upper() if match else None


def suggested_role(series_description: Any, modality: Any) -> str:
    description = str(series_description or "").strip().lower()
    modality = str(modality or "").strip().upper()
    if "analysis mask" in description:
        return "FTV_MASK_REVIEW"
    if "volser" in description and re.search(r"\bpe\d+\b", description):
        return "PE_MAP_REVIEW"
    if "volser" in description and re.search(r"\bser\b", description):
        return "SER_MAP_REVIEW"
    if "volser" in description and "original dce" in description:
        return "CROPPED_DCE_ORIGINAL_REVIEW"
    if modality == "MR" and ("t2" in description or "t2" in str(series_description or "").lower()):
        return "T2_ANATOMICAL_CANDIDATE"
    if modality == "SEG":
        return "SEG_REVIEW"
    return "MANUAL_SERIES_REVIEW"


def build_candidate_rows(
    series_rows: list[dict[str, Any]],
    study_tags: dict[str, dict[str, Any]],
    target_subject_ids: set[str],
    collection: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    records: list[dict[str, Any]] = []
    matched_studies: set[str] = set()
    matched_subjects: set[str] = set()
    blocked_studies = 0
    for row in series_rows:
        study_uid = str(row.get("StudyInstanceUID") or row.get("studyInstanceUID") or "").strip()
        series_uid = str(row.get("SeriesInstanceUID") or row.get("seriesInstanceUID") or "").strip()
        if not study_uid or not series_uid:
            continue
        tags = study_tags.get(study_uid, {})
        raw_subject = tags.get("trial_subject_id")
        subject_id = canonical_subject_id(raw_subject)
        if (
            subject_id is None
            or tags.get("subject_tag_value_count") != 1
            or tags.get("visit_tag_value_count") != 1
        ):
            blocked_studies += 1
            continue
        if subject_id not in target_subject_ids:
            continue
        visit_code = str(tags.get("visit_code") or "")
        modality = row.get("Modality") or row.get("modality")
        description = row.get("SeriesDescription") or row.get("seriesDescription")
        records.append({
            "trial_subject_id": subject_id,
            "visit_code": visit_code,
            "visit_stage": canonical_stage(visit_code),
            "collection": collection,
            "patient_id": row.get("PatientID") or row.get("patientID"),
            "study_instance_uid": study_uid,
            "series_instance_uid": series_uid,
            "modality": modality,
            "series_number": row.get("SeriesNumber") or row.get("seriesNumber"),
            "series_description": description,
            "protocol_name": row.get("ProtocolName") or row.get("protocolName"),
            "manufacturer": row.get("Manufacturer") or row.get("manufacturer"),
            "manufacturer_model_name": row.get("ManufacturerModelName") or row.get("manufacturerModelName"),
            "image_count": row.get("ImageCount") or row.get("imageCount"),
            "suggested_role": suggested_role(description, modality),
            "include": 0,
            "review_note": "Candidate only; confirm visit, anatomy and series content before inclusion",
        })
        matched_studies.add(study_uid)
        matched_subjects.add(subject_id)

    columns = [
        "trial_subject_id", "visit_code", "visit_stage", "collection", "patient_id",
        "study_instance_uid", "series_instance_uid", "modality", "series_number",
        "series_description", "protocol_name", "manufacturer", "manufacturer_model_name",
        "image_count", "suggested_role", "include", "review_note",
    ]
    frame = pd.DataFrame(records, columns=columns)
    summary = {
        "collection": collection,
        "target_subject_count": len(target_subject_ids),
        "matched_subject_count": len(matched_subjects),
        "unmatched_target_subject_count": len(target_subject_ids - matched_subjects),
        "matched_study_count": len(matched_studies),
        "candidate_series_count": len(frame),
        "candidate_visit_stage_counts": (
            frame["visit_stage"].fillna("UNKNOWN").value_counts().sort_index().astype(int).to_dict()
            if not frame.empty else {}
        ),
        "suggested_role_counts": (
            frame["suggested_role"].value_counts().sort_index().astype(int).to_dict()
            if not frame.empty else {}
        ),
        "blocked_study_count_missing_or_ambiguous_trial_tags": blocked_studies,
        "subject_mapping": "Verified from DICOM ClinicalTrialSubjectID (0012,0040), not inferred from PatientID",
        "selection_policy": "include defaults to 0; human review is mandatory",
        "privacy": "Candidate CSV contains participant IDs and SeriesInstanceUIDs; keep local and private",
    }
    return frame, summary


def query_subset(
    feature_workbook: Path,
    sheet_name: str,
    id_column: str,
    collection: str,
    workers: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    feature = pd.read_excel(feature_workbook, sheet_name=sheet_name)
    if id_column not in feature.columns:
        raise ValueError(f"Feature workbook is missing identifier column {id_column!r}")
    raw_ids = feature[id_column].dropna().tolist()
    target_ids = {canonical_subject_id(value) for value in raw_ids}
    if None in target_ids:
        raise ValueError("All feature workbook IDs must resolve to six-digit trial subject IDs")
    target_ids = {value for value in target_ids if value is not None}
    if len(target_ids) != len(raw_ids):
        raise ValueError("Feature workbook contains duplicate or non-unique subject IDs")

    series_rows = fetch_public_series(collection)
    study_to_series: dict[str, list[dict[str, Any]]] = {}
    for row in series_rows:
        study_uid = str(row.get("StudyInstanceUID") or row.get("studyInstanceUID") or "").strip()
        series_uid = str(row.get("SeriesInstanceUID") or row.get("seriesInstanceUID") or "").strip()
        if study_uid and series_uid:
            study_to_series.setdefault(study_uid, []).append(row)

    token = request_guest_token()
    representative: dict[str, str] = {}
    for study_uid, rows in study_to_series.items():
        rows = sorted(
            rows,
            key=lambda r: (
                str(r.get("Modality") or r.get("modality") or "") != "MR",
                str(r.get("SeriesNumber") or r.get("seriesNumber") or ""),
                str(r.get("SeriesInstanceUID") or r.get("seriesInstanceUID") or ""),
            ),
        )
        representative[study_uid] = str(
            rows[0].get("SeriesInstanceUID") or rows[0].get("seriesInstanceUID")
        )

    study_tags: dict[str, dict[str, Any]] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        futures = {
            pool.submit(fetch_study_trial_tags, uid, token): study_uid
            for study_uid, uid in representative.items()
        }
        for future in concurrent.futures.as_completed(futures):
            study_uid = futures[future]
            study_tags[study_uid] = future.result()

    frame, summary = build_candidate_rows(series_rows, study_tags, target_ids, collection)
    summary["api_series_rows_returned"] = len(series_rows)
    summary["studies_with_trial_tags_queried"] = len(study_tags)
    summary["study_tag_query_workers"] = max(1, workers)
    return frame, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature-workbook", required=True)
    parser.add_argument("--sheet-name", default="datawith4visits")
    parser.add_argument("--id-column", default="CLINICAL-TRIAL-SUBJECT-ID")
    parser.add_argument("--collection", default="ISPY2")
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--candidate-csv", required=True)
    parser.add_argument("--summary-json", required=True)
    args = parser.parse_args()
    if args.workers < 1 or args.workers > 6:
        parser.error("--workers must be between 1 and 6 to avoid aggressive API traffic")
    frame, summary = query_subset(
        Path(args.feature_workbook), args.sheet_name, args.id_column,
        args.collection, args.workers,
    )
    candidate_path = Path(args.candidate_csv)
    summary_path = Path(args.summary_json)
    candidate_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(candidate_path, index=False)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
