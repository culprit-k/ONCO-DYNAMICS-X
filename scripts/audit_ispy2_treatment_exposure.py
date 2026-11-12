"""Align administered drug events to TCIA-shifted MRI visit dates.

The script refuses original-calendar treatment dates. Dates must already be
expressed on the same TCIA modified longitudinal date axis as the DICOM scans.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import pandas as pd

REQUIRED_EVENT_COLUMNS = {
    "patient_id", "agent", "administration_datetime", "dose", "dose_unit",
    "route", "administration_status", "time_basis", "source_reference",
}
REQUIRED_VISIT_COLUMNS = {
    "patient_id", "visit_code", "scan_datetime", "longitudinal_date_flag",
}


def align_exposures(events: pd.DataFrame, visits: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    missing_events = REQUIRED_EVENT_COLUMNS - set(events.columns)
    missing_visits = REQUIRED_VISIT_COLUMNS - set(visits.columns)
    if missing_events or missing_visits:
        raise ValueError(
            f"Missing columns; events={sorted(missing_events)}, visits={sorted(missing_visits)}"
        )
    e = events.copy()
    v = visits.copy()
    for frame in (e, v):
        frame["patient_id"] = frame["patient_id"].astype("string").str.strip().str.zfill(6)
    if e["patient_id"].isna().any() or v["patient_id"].isna().any():
        raise ValueError("Patient IDs must be present")
    if not e["time_basis"].astype(str).str.upper().eq("TCIA_MODIFIED_LONGITUDINAL").all():
        raise ValueError(
            "Treatment dates are not verified on the TCIA modified longitudinal date axis. "
            "Original-calendar dates cannot be aligned without a documented patient-specific shift."
        )
    flags = set(v["longitudinal_date_flag"].fillna("").astype(str).str.upper())
    if not flags or not flags.issubset({"MODIFIED", "UNMODIFIED"}):
        raise ValueError("Visit DICOM date-modification flag is missing or inconsistent")
    if "MODIFIED" not in flags:
        raise ValueError("No TCIA MODIFIED DICOM dates detected in visit timeline")

    e["administration_status"] = e["administration_status"].astype(str).str.upper().str.strip()
    excluded_status_counts = dict(Counter(e.loc[
        e["administration_status"] != "ADMINISTERED", "administration_status"
    ]))
    e = e[e["administration_status"] == "ADMINISTERED"].copy()
    if e.empty:
        raise ValueError("No administered treatment events were supplied")
    if e[["agent", "dose_unit", "route", "source_reference"]].isna().any().any():
        raise ValueError("Actual administration records require agent, dose unit, route and source")
    e["agent"] = e["agent"].astype(str).str.strip()
    e["dose_unit"] = e["dose_unit"].astype(str).str.strip()
    e["dose"] = pd.to_numeric(e["dose"], errors="raise")
    if (~e["dose"].map(lambda x: pd.notna(x) and x >= 0)).any():
        raise ValueError("Dose must be finite and non-negative")
    e["administration_datetime"] = pd.to_datetime(
        e["administration_datetime"], errors="raise", utc=True
    )
    v["scan_datetime"] = pd.to_datetime(v["scan_datetime"], errors="raise", utc=True)
    if e["administration_datetime"].isna().any() or v["scan_datetime"].isna().any():
        raise ValueError("Administration and scan dates must be parseable")

    output = []
    matched_visits = 0
    unmatched_visits = 0
    for row in v.itertuples(index=False):
        patient_events = e[e["patient_id"] == row.patient_id]
        if patient_events.empty:
            unmatched_visits += 1
            continue
        prior = patient_events[
            patient_events["administration_datetime"] <= row.scan_datetime
        ]
        if prior.empty:
            unmatched_visits += 1
            continue
        for agent, agent_events in prior.groupby("agent"):
            units = agent_events["dose_unit"].unique()
            if len(units) != 1:
                raise ValueError(
                    f"Mixed dose units for patient {row.patient_id}, agent {agent}; "
                    "harmonise units before cumulative-dose calculation"
                )
            agent_events = agent_events.sort_values("administration_datetime")
            last = agent_events.iloc[-1]
            days_since = (
                row.scan_datetime - last["administration_datetime"]
            ).total_seconds() / 86400.0
            output.append({
                "patient_id": row.patient_id,
                "visit_code": row.visit_code,
                "scan_datetime": row.scan_datetime.isoformat(),
                "agent": agent,
                "last_administration_datetime": last["administration_datetime"].isoformat(),
                "days_since_last_administration": float(days_since),
                "cumulative_administered_dose_before_scan": float(agent_events["dose"].sum()),
                "dose_unit": units[0],
                "last_route": str(last["route"]),
                "administered_event_count_before_scan": int(len(agent_events)),
                "source_time_basis": "TCIA_MODIFIED_LONGITUDINAL",
            })
        matched_visits += 1

    aligned = pd.DataFrame(output)
    summary = {
        "visit_rows": int(len(v)),
        "visits_with_at_least_one_prior_administration": int(matched_visits),
        "visits_without_prior_administration": int(unmatched_visits),
        "administered_event_rows_used": int(len(e)),
        "patients_with_events": int(e["patient_id"].nunique()),
        "agents": sorted(e["agent"].unique().tolist()),
        "excluded_non_administered_status_counts": excluded_status_counts,
        "time_basis": "TCIA_MODIFIED_LONGITUDINAL",
        "interpretation": (
            "Days-since-administration and cumulative dose are descriptive alignment features. "
            "They are not pharmacokinetic concentrations or mechanistic drug-effect parameters."
        ),
        "privacy": "Aligned CSV contains patient IDs, shifted scan dates and treatment event data; keep private.",
    }
    return aligned, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events-csv", required=True, help="Private actual administration event table")
    parser.add_argument("--visits-csv", required=True, help="Private DICOM-derived visit timeline")
    parser.add_argument("--aligned-csv", required=True, help="Private local aligned event table")
    parser.add_argument("--summary-json", required=True, help="Aggregate-only summary JSON")
    args = parser.parse_args()
    events = pd.read_csv(args.events_csv, dtype={"patient_id": "string"})
    visits = pd.read_csv(args.visits_csv, dtype={"patient_id": "string"})
    aligned, summary = align_exposures(events, visits)
    output = Path(args.aligned_csv)
    summary_path = Path(args.summary_json)
    output.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    aligned.to_csv(output, index=False)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
