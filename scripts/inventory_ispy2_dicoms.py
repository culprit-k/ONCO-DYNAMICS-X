"""Create a local-only, metadata-only DICOM inventory for I-SPY 2.

No pixel data are read and no patient-level inventory is intended for GitHub.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd

from onco_dynamics_x.imaging.dicom_metadata import extract_series_metadata

def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, default=str)


def inventory_dicom_tree(root: Path, *, max_files: int | None = None) -> tuple[pd.DataFrame, dict]:
    try:
        import pydicom
    except ImportError as exc:  # pragma: no cover - exercised in environment setup
        raise RuntimeError("Install the optional imaging dependencies: pip install -e '.[imaging]'") from exc

    if not root.exists() or not root.is_dir():
        raise ValueError(f"DICOM root is not a directory: {root}")

    groups: dict[tuple, dict] = {}
    inspected = 0
    dicom_count = 0
    errors = 0
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        if max_files is not None and inspected >= max_files:
            break
        inspected += 1
        try:
            ds = pydicom.dcmread(str(path), stop_before_pixels=True, force=True)
        except Exception:
            errors += 1
            continue
        # Force parsing can open non-DICOM files. Require essential DICOM identity tags.
        metadata = extract_series_metadata(ds)
        series_uid = metadata.get("series_instance_uid")
        study_uid = metadata.get("study_instance_uid")
        if not series_uid or not study_uid or not metadata.get("modality"):
            continue
        dicom_count += 1
        key = (
            metadata.get("patient_id") or "",
            metadata.get("visit_code") or "",
            study_uid,
            series_uid,
        )
        record = groups.get(key)
        if record is None:
            record = dict(metadata)
            record.update({
                "series_directory": str(path.parent.resolve()),
                "file_paths": [],
                "file_count": 0,
                "sop_uids": set(),
                "study_datetimes": set(),
                "acquisition_datetimes": [],
                "date_flags": set(),
                "patient_id_sources": set(),
                "ftv_results_by_signature": {},
                "ser_parameters_by_signature": {},
            })
            groups[key] = record

        record["file_count"] += 1
        record["file_paths"].append(str(path.resolve()))
        if metadata.get("sop_instance_uid"):
            record["sop_uids"].add(metadata["sop_instance_uid"])
        if metadata.get("study_datetime"):
            record["study_datetimes"].add(metadata["study_datetime"])
        if metadata.get("acquisition_datetime"):
            record["acquisition_datetimes"].append(metadata["acquisition_datetime"])
        if metadata.get("longitudinal_date_flag"):
            record["date_flags"].add(metadata["longitudinal_date_flag"].upper())
        if metadata.get("patient_id_source"):
            record["patient_id_sources"].add(metadata["patient_id_source"])
        if metadata.get("ftv_results"):
            signature = _json(metadata["ftv_results"])
            record["ftv_results_by_signature"][signature] = metadata["ftv_results"]
        if metadata.get("ser_analysis_parameters"):
            signature = _json(metadata["ser_analysis_parameters"])
            record["ser_parameters_by_signature"][signature] = metadata["ser_analysis_parameters"]

    rows = []
    for record in groups.values():
        result = dict(record)
        acq = sorted(datetime.fromisoformat(x) for x in record["acquisition_datetimes"])
        study_dates = sorted(record["study_datetimes"])
        flags = sorted(record["date_flags"])
        result["study_datetime"] = study_dates[0] if len(study_dates) == 1 else None
        result["study_datetime_consistent"] = len(study_dates) == 1
        result["scan_datetime"] = result["study_datetime"] or (
            acq[0].isoformat() if acq else None
        )
        result["scan_datetime_source"] = (
            "StudyDate/StudyTime" if result["study_datetime"] else
            (record.get("acquisition_datetime_source") if acq else None)
        )
        result["acquisition_start_datetime"] = acq[0].isoformat() if acq else None
        result["acquisition_end_datetime"] = acq[-1].isoformat() if acq else None
        result["acquisition_datetime_count"] = len(acq)
        result["date_flags"] = "|".join(flags)
        result["date_flag_consistent"] = len(flags) == 1
        result["patient_id_sources"] = "|".join(sorted(record["patient_id_sources"]))
        result["unique_sop_instance_count"] = len(record["sop_uids"])
        result["study_datetime_values_json"] = _json(study_dates)
        result["file_paths_json"] = _json(sorted(set(record["file_paths"])))
        result["ftv_results_json"] = _json(list(record["ftv_results_by_signature"].values()))
        result["ser_analysis_parameters_json"] = _json(list(record["ser_parameters_by_signature"].values()))
        for key in (
            "sop_uids", "study_datetimes", "acquisition_datetimes", "date_flags",
            "patient_id_sources", "ftv_results_by_signature", "ser_parameters_by_signature",
            "file_paths",
        ):
            result.pop(key, None)
        rows.append(result)

    df = pd.DataFrame(rows)
    visit_counts = Counter(
        row.get("visit_code") or "MISSING" for row in rows
    )
    summary = {
        "input_root_name": root.name,
        "files_examined": inspected,
        "dicom_image_or_segmentation_files": dicom_count,
        "series_count": len(rows),
        "unique_patient_count": int(df["patient_id"].dropna().nunique()) if not df.empty else 0,
        "series_without_patient_id": int(df["patient_id"].isna().sum()) if not df.empty else 0,
        "series_without_visit_code": int((df["visit_code"].isna() | (df["visit_code"] == "")).sum()) if not df.empty else 0,
        "visit_code_series_counts": dict(visit_counts),
        "series_with_scan_datetime": int(df["scan_datetime"].notna().sum()) if not df.empty else 0,
        "series_with_longitudinal_date_flag": int(df["date_flags"].astype(bool).sum()) if not df.empty else 0,
        "series_with_ftv_private_results": int(
            df["ftv_results_json"].fillna("[]").ne("[]").sum()
        ) if not df.empty else 0,
        "series_with_ser_analysis_parameters": int(
            df["ser_analysis_parameters_json"].fillna("[]").ne("[]").sum()
        ) if not df.empty else 0,
        "parse_errors": errors,
        "privacy": "CSV contains patient IDs, DICOM UIDs, shifted dates and local paths; keep private.",
        "date_warning": (
            "TCIA DICOM date fields can be modified. Inspect LongitudinalTemporalInformationModified "
            "(0028,0303); never treat shifted dates as real calendar dates."
        ),
        "timing_warning": (
            "DCE private timing fields are core-lab estimates; the collection guide states that "
            "contrast-injection timing accuracy cannot be guaranteed."
        ),
    }
    return df, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dicom-root", required=True, help="Local TCIA DICOM directory")
    parser.add_argument("--inventory-csv", required=True, help="Private local series inventory CSV")
    parser.add_argument("--summary-json", required=True, help="Aggregate-only summary JSON")
    parser.add_argument("--max-files", type=int, help="Optional small test cap; not for production audit")
    args = parser.parse_args()
    if args.max_files is not None and args.max_files < 1:
        parser.error("--max-files must be positive")
    df, summary = inventory_dicom_tree(Path(args.dicom_root), max_files=args.max_files)
    inventory_path = Path(args.inventory_csv)
    summary_path = Path(args.summary_json)
    inventory_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(inventory_path, index=False)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
