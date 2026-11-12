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


def select_inventory_rows(inventory: pd.DataFrame, selection: pd.DataFrame) -> pd.DataFrame:
    """Filter inventory to exact human-reviewed patient/stage/series triples."""
    required_inventory = {"patient_id", "series_instance_uid"}
    required_selection = {"patient_id", "visit_stage", "series_instance_uid", "include"}
    missing_inventory = required_inventory - set(inventory.columns)
    missing_selection = required_selection - set(selection.columns)
    if missing_inventory or missing_selection:
        raise ValueError(
            f"Missing selection columns; inventory={sorted(missing_inventory)}, "
            f"selection={sorted(missing_selection)}"
        )
    data = inventory.copy()
    chosen = selection.copy()
    for frame in (data, chosen):
        frame["patient_id"] = frame["patient_id"].astype("string").str.strip().str.zfill(6)
        frame["series_instance_uid"] = frame["series_instance_uid"].astype("string").str.strip()
    if "visit_stage" not in data.columns:
        data["visit_stage"] = data["visit_code"].astype("string").str.upper().str.extract(
            r"^(T[0-3])(?:[-_].*)?$", expand=False
        )
    data["visit_stage"] = data["visit_stage"].astype("string").str.upper()
    chosen["visit_stage"] = chosen["visit_stage"].astype("string").str.upper()
    raw_include = chosen["include"].astype(str).str.strip().str.lower()
    valid = {"1", "0", "true", "false", "yes", "no"}
    if not set(raw_include.unique()).issubset(valid):
        raise ValueError("Selection include values must be 1/0, true/false or yes/no")
    chosen = chosen.loc[raw_include.isin({"1", "true", "yes"})].copy()
    if chosen.empty:
        raise ValueError("No patient/visit/series rows were selected")
    if chosen["series_instance_uid"].duplicated().any():
        raise ValueError("Selected series contains duplicate SeriesInstanceUID values")
    triples = ["patient_id", "visit_stage", "series_instance_uid"]
    check = chosen[triples].merge(
        data[triples].drop_duplicates(),
        on=triples,
        how="left",
        validate="one_to_one",
        indicator=True,
    )
    if not check["_merge"].eq("both").all():
        raise ValueError("At least one selected patient/visit/series triple is absent from inventory")
    return data[data["series_instance_uid"].isin(set(chosen["series_instance_uid"]))].copy()


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
    if "visit_stage" in data.columns:
        data["_visit_stage"] = data["visit_stage"].astype("string").str.upper()
    else:
        data["_visit_stage"] = data["visit_code"].str.extract(r"^(T[0-3])(?:[-_].*)?$", expand=False)
    data = data[data["_visit_stage"].isin(VISIT_ORDER)].copy()
    if data.empty:
        raise ValueError("No T0/T1/T2/T3 DICOM visit records found")

    records = []
    for (patient_id, visit, study_uid), group in data.groupby(
        ["patient_id", "_visit_stage", "study_instance_uid"], dropna=False
    ):
        raw_codes = sorted(set(group["visit_code"].dropna().astype(str)))
        study_datetimes = (
            sorted(set(group["study_datetime"].dropna().astype(str)))
            if "study_datetime" in group.columns else []
        )
        inconsistent_study_rows = False
        if "study_datetime_consistent" in group.columns:
            for _, series_row in group.iterrows():
                has_study_values = False
                if "study_datetime_values_json" in group.columns:
                    raw_values = series_row.get("study_datetime_values_json")
                    has_study_values = (
                        isinstance(raw_values, str)
                        and raw_values not in ("", "[]", "null")
                    )
                if has_study_values and not bool(series_row["study_datetime_consistent"]):
                    inconsistent_study_rows = True

        if len(study_datetimes) == 1 and not inconsistent_study_rows:
            scan_datetime = study_datetimes[0]
            scan_datetime_source = "StudyDate/StudyTime"
            precision_values = (
                sorted(set(group["study_time_precision"].dropna().astype(str)))
                if "study_time_precision" in group.columns else ["date_and_time"]
            )
        elif not study_datetimes and not inconsistent_study_rows:
            acquisition_starts = (
                sorted(set(group["acquisition_start_datetime"].dropna().astype(str)))
                if "acquisition_start_datetime" in group.columns else []
            )
            if not acquisition_starts:
                fallback_dates = sorted(set(group["scan_datetime"].dropna().astype(str)))
                acquisition_starts = fallback_dates
            scan_datetime = acquisition_starts[0] if acquisition_starts else None
            scan_datetime_source = "earliest available series acquisition datetime"
            precision_values = (
                sorted(set(group["acquisition_time_precision"].dropna().astype(str)))
                if "acquisition_time_precision" in group.columns else ["date_and_time"]
            )
        else:
            scan_datetime = None
            scan_datetime_source = "ambiguous study/acquisition date metadata"
            precision_values = []

        date_flags = sorted({
            token
            for raw in group["date_flags"].fillna("").astype(str)
            for token in raw.split("|")
            if token
        })
        date_flag_valid = (
            bool(group["date_flag_consistent"].fillna(False).astype(bool).all())
            if "date_flag_consistent" in group.columns
            else len(date_flags) == 1 and date_flags[0] in VALID_DATE_FLAGS
        )
        records.append({
            "patient_id": str(patient_id),
            "visit_code": str(visit),
            "visit_code_raw": "|".join(raw_codes),
            "study_instance_uid": "" if pd.isna(study_uid) else str(study_uid),
            "scan_datetime": scan_datetime,
            "scan_datetime_consistent": scan_datetime is not None,
            "scan_datetime_source": scan_datetime_source,
            "scan_time_precision": "|".join(precision_values),
            "longitudinal_date_flag": "|".join(date_flags),
            "date_flag_valid": date_flag_valid,
            "series_count": int(len(group)),
            "study_datetime_consistent_across_series": not inconsistent_study_rows,
        })

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
            continue
        t0 = patient[patient["visit_code"] == "T0"]
        if len(t0) != 1:
            visits.loc[idx, "interval_status"] = "blocked_missing_or_multiple_T0_studies"
            continue
        if not bool(t0.iloc[0]["scan_datetime_consistent"]):
            visits.loc[idx, "interval_status"] = "blocked_inconsistent_T0_datetime"
            continue
        t0_datetime = str(t0.iloc[0]["scan_datetime"])
        try:
            # TCIA-shifted dates are used only for within-patient elapsed intervals.
            pd.Timestamp(t0_datetime)
        except Exception:
            visits.loc[idx, "interval_status"] = "blocked_unparseable_T0_datetime"
            continue
        t0_precision = str(t0.iloc[0]["scan_time_precision"])
        if "|" in t0_precision:
            visits.loc[idx, "interval_status"] = "blocked_mixed_T0_time_precision"
            continue

        ordered = patient.copy()
        ordered["_rank"] = ordered["visit_code"].map(VISIT_ORDER)
        ordered = ordered.sort_values("_rank")
        previous_valid_datetime = None
        previous_rank = None
        for row_index, row in ordered.iterrows():
            visit = str(row["visit_code"])
            if pd.isna(row["scan_datetime"]):
                visits.loc[row_index, "interval_status"] = "blocked_missing_visit_datetime"
                previous_valid_datetime = None
                previous_rank = None
                continue
            if not bool(row["scan_datetime_consistent"]):
                visits.loc[row_index, "interval_status"] = "blocked_inconsistent_visit_datetime"
                previous_valid_datetime = None
                previous_rank = None
                continue
            if not bool(row["study_datetime_consistent_across_series"]):
                visits.loc[row_index, "interval_status"] = "blocked_inconsistent_study_datetime_metadata"
                previous_valid_datetime = None
                previous_rank = None
                continue
            current = str(row["scan_datetime"])
            current_precision = str(row["scan_time_precision"])
            if "|" in current_precision or current_precision not in {"date_only", "date_and_time"}:
                visits.loc[row_index, "interval_status"] = "blocked_unknown_time_precision"
                previous_valid_datetime = None
                previous_rank = None
                continue
            if current_precision == "date_and_time" and t0_precision == "date_and_time":
                days_from_t0 = elapsed_days(t0_datetime, current)
                interval_precision = "date_and_time"
            else:
                days_from_t0 = float(
                    (pd.Timestamp(current).date() - pd.Timestamp(t0_datetime).date()).days
                )
                interval_precision = "date"
            if days_from_t0 < 0:
                visits.loc[row_index, "interval_status"] = "blocked_negative_time_from_T0"
                previous_valid_datetime = None
                previous_rank = None
                continue
            visits.loc[row_index, "elapsed_days_from_T0"] = days_from_t0
            visits.loc[row_index, "elapsed_time_precision"] = interval_precision
            if visit == "T0":
                visits.loc[row_index, "days_from_previous_stage"] = 0.0
                visits.loc[row_index, "interval_status"] = "baseline"
            elif previous_valid_datetime is not None and previous_rank is not None:
                if VISIT_ORDER[visit] > previous_rank:
                    previous_row = ordered.loc[ordered["_rank"] == previous_rank].iloc[-1]
                    previous_precision = str(previous_row["scan_time_precision"])
                    if previous_precision == "date_and_time" and current_precision == "date_and_time":
                        stage_interval = elapsed_days(previous_valid_datetime, current)
                        stage_precision = "date_and_time"
                    else:
                        stage_interval = float(
                            (pd.Timestamp(current).date() - pd.Timestamp(previous_valid_datetime).date()).days
                        )
                        stage_precision = "date"
                    visits.loc[row_index, "days_from_previous_stage"] = stage_interval
                    visits.loc[row_index, "stage_interval_precision"] = stage_precision
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
        "elapsed_time_precision_counts": dict(Counter(visits.get("elapsed_time_precision", pd.Series(dtype=str)).dropna())) ,
        "stage_interval_precision_counts": dict(Counter(visits.get("stage_interval_precision", pd.Series(dtype=str)).dropna())),
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
    parser.add_argument(
        "--selection-csv",
        help="Optional reviewed patient/visit/series selection CSV; resolves test/retest ambiguity",
    )
    parser.add_argument("--timeline-csv", required=True, help="Private local patient-visit timeline")
    parser.add_argument("--summary-json", required=True, help="Aggregate-only summary")
    args = parser.parse_args()
    inventory = pd.read_csv(args.inventory_csv, dtype={"patient_id": "string"})
    if args.selection_csv:
        selection = pd.read_csv(args.selection_csv, dtype={"patient_id": "string"})
        inventory = select_inventory_rows(inventory, selection)
    timeline, summary = build_visit_timeline(inventory)
    summary["series_selection_applied"] = bool(args.selection_csv)
    summary["series_selection_policy"] = (
        "exact reviewed patient/visit/series triples"
        if args.selection_csv else "none; ambiguous repeated studies are blocked"
    )
    output = Path(args.timeline_csv)
    summary_path = Path(args.summary_json)
    output.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    timeline.to_csv(output, index=False)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
