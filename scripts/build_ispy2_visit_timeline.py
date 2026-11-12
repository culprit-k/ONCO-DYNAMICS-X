"""Build patient-level I-SPY 2 visit intervals from TCIA-modified DICOM dates.

Only relative day differences are calculated when the DICOM date-modification
flag is present and consistent. Absolute calendar dates are never interpreted
as real-world dates.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import pandas as pd

from onco_dynamics_x.imaging.dicom_metadata import elapsed_days

VISIT_ORDER = {"T0": 0, "T1": 1, "T2": 2, "T3": 3}
VALID_DATE_FLAGS = {"MODIFIED", "UNMODIFIED"}


def build_visit_timeline(inventory: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    required = {
        "patient_id", "visit_code", "study_instance_uid", "scan_datetime",
        "scan_datetime_source", "date_flags", "study_datetime_consistent",
    }
    missing = required - set(inventory.columns)
    if missing:
        raise ValueError(f"Inventory missing required columns: {sorted(missing)}")

    data = inventory.copy()
    data["patient_id"] = data["patient_id"].astype("string")
    data["visit_code"] = data["visit_code"].astype("string").str.upper()
    data = data[data["visit_code"].isin(VISIT_ORDER)].copy()
    if data.empty:
        raise ValueError("No T0/T1/T2/T3 DICOM visit records found")

    records = []
    for (patient_id, visit, study_uid), group in data.groupby(
        ["patient_id", "visit_code", "study_instance_uid"], dropna=False
    ):
        dates = sorted(set(group["scan_datetime"].dropna().astype(str)))
        flag_values = set()
        for raw in group["date_flags"].fillna("").astype(str):
            flag_values.update(x for x in raw.split("|") if x)
        date_flags = sorted(flag_values)
        record = {
            "patient_id": str(patient_id),
            "visit_code": str(visit),
            "study_instance_uid": "" if pd.isna(study_uid) else str(study_uid),
            "scan_datetime": dates[0] if len(dates) == 1 else None,
            "scan_datetime_consistent": len(dates) == 1,
            "scan_datetime_source": "|".join(sorted(set(
                group["scan_datetime_source"].dropna().astype(str)
            ))),
            "longitudinal_date_flag": "|".join(date_flags),
            "date_flag_valid": len(date_flags) == 1 and date_flags[0] in VALID_DATE_FLAGS,
            "series_count": int(len(group)),
            "study_datetime_consistent_across_series": bool(
                group["study_datetime_consistent"].fillna(False).astype(bool).all()
            ),
        }
        records.append(record)

    visits = pd.DataFrame(records)
    visits["elapsed_days_from_T0"] = pd.NA
    visits["days_from_previous_stage"] = pd.NA
    visits["interval_status"] = "not_evaluated"

    for patient_id, indices in visits.groupby("patient_id").groups.items():
        idx = list(indices)
        patient = visits.loc[idx].copy()
        flag_set = set(
            value
            for raw in patient["longitudinal_date_flag"].astype(str)
            for value in raw.split("|")
            if value
        )
        all_dates_valid = (
            len(flag_set) == 1
            and next(iter(flag_set)) in VALID_DATE_FLAGS
            and patient["date_flag_valid"].all()
        )
        if not all_dates_valid:
            visits.loc[idx, "interval_status"] = "blocked_missing_or_inconsistent_date_flag"
            continue

        ambiguous = patient.groupby("visit_code").size()
        if (ambiguous > 1).any():
            visits.loc[idx, "interval_status"] = "blocked_multiple_studies_within_visit"
        t0 = patient[patient["visit_code"] == "T0"]
        if len(t0) != 1:
            visits.loc[idx, "interval_status"] = "blocked_missing_or_multiple_T0_studies"
            continue
        if not bool(t0.iloc[0]["scan_datetime_consistent"]):
            visits.loc[idx, "interval_status"] = "blocked_inconsistent_T0_datetime"
            continue
        t0_datetime = str(t0.iloc[0]["scan_datetime"])
        try:
            # Validate parseability without converting to a real calendar date.
            pd.Timestamp(t0_datetime)
        except Exception:
            visits.loc[idx, "interval_status"] = "blocked_unparseable_T0_datetime"
            continue

        ordered = patient.copy()
        ordered["_rank"] = ordered["visit_code"].map(VISIT_ORDER)
        ordered = ordered.sort_values("_rank")
        previous_valid_datetime = None
        previous_rank = None
        for row_index, row in ordered.iterrows():
            visit = str(row["visit_code"])
            if not bool(row["scan_datetime_consistent"]) or not row["scan_datetime"]:
                visits.loc[row_index, "interval_status"] = "blocked_inconsistent_visit_datetime"
                previous_valid_datetime = None
                previous_rank = None
                continue
            current = str(row["scan_datetime"])
            days_from_t0 = elapsed_days(t0_datetime, current)
            if days_from_t0 < 0:
                visits.loc[row_index, "interval_status"] = "blocked_negative_time_from_T0"
                continue
            visits.loc[row_index, "elapsed_days_from_T0"] = days_from_t0
            if visit == "T0":
                visits.loc[row_index, "days_from_previous_stage"] = 0.0
                visits.loc[row_index, "interval_status"] = "baseline"
            elif previous_valid_datetime is not None and previous_rank is not None:
                if VISIT_ORDER[visit] > previous_rank:
                    visits.loc[row_index, "days_from_previous_stage"] = elapsed_days(
                        previous_valid_datetime, current
                    )
                    visits.loc[row_index, "interval_status"] = "relative_interval_valid"
                else:
                    visits.loc[row_index, "interval_status"] = "blocked_nonmonotone_visit_order"
            else:
                visits.loc[row_index, "interval_status"] = "T0_relative_days_valid_previous_stage_ambiguous"
            previous_valid_datetime = current
            previous_rank = VISIT_ORDER[visit]

    summary = {
        "patient_visit_studies": int(len(visits)),
        "unique_patients": int(visits["patient_id"].nunique()),
        "visit_study_counts": dict(Counter(visits["visit_code"])),
        "studies_with_scan_datetime": int(visits["scan_datetime"].notna().sum()),
        "studies_with_valid_date_flag": int(visits["date_flag_valid"].sum()),
        "visits_with_t0_relative_days": int(visits["elapsed_days_from_T0"].notna().sum()),
        "visits_with_stage_intervals": int(visits["days_from_previous_stage"].notna().sum()),
        "interval_status_counts": dict(Counter(visits["interval_status"])),
        "date_semantics": (
            "DICOM calendar dates are privacy-shifted where the MODIFIED flag is present; "
            "only within-patient elapsed intervals are used."
        ),
        "privacy": "Output CSV contains participant IDs, StudyInstanceUIDs and shifted dates; keep private.",
    }
    return visits, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory-csv", required=True)
    parser.add_argument("--timeline-csv", required=True, help="Private local patient-visit timeline")
    parser.add_argument("--summary-json", required=True, help="Aggregate-only summary")
    args = parser.parse_args()
    inventory = pd.read_csv(args.inventory_csv, dtype={"patient_id": "string"})
    timeline, summary = build_visit_timeline(inventory)
    output = Path(args.timeline_csv)
    summary_path = Path(args.summary_json)
    output.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    timeline.to_csv(output, index=False)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
