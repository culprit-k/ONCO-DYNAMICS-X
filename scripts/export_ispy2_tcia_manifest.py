"""Export a reviewed TCIA SeriesInstanceUID manifest from private candidate rows."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

REQUIRED_COLUMNS = {
    "trial_subject_id",
    "visit_code",
    "visit_stage",
    "series_instance_uid",
    "modality",
    "series_description",
    "include",
    "review_note",
}
DEFAULT_REVIEW_NOTE = (
    "Candidate only; confirm visit, anatomy and series content before inclusion"
)


def export_reviewed_manifest(candidates: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    missing = REQUIRED_COLUMNS - set(candidates.columns)
    if missing:
        raise ValueError(f"Candidate file is missing columns: {sorted(missing)}")
    frame = candidates.copy()
    flags = frame["include"].astype("string").str.strip().str.lower()
    valid = {"0", "1", "true", "false", "yes", "no"}
    if not set(flags.dropna().unique()).issubset(valid) or flags.isna().any():
        raise ValueError("include must use 1/0, true/false or yes/no")
    selected = frame.loc[flags.isin({"1", "true", "yes"})].copy()
    if selected.empty:
        raise ValueError("No series marked include=1")
    if selected["series_instance_uid"].isna().any():
        raise ValueError("Selected series must have SeriesInstanceUID")
    selected["series_instance_uid"] = (
        selected["series_instance_uid"].astype("string").str.strip()
    )
    if selected["series_instance_uid"].eq("").any():
        raise ValueError("Selected series must have non-empty SeriesInstanceUID")
    if selected["series_instance_uid"].duplicated().any():
        raise ValueError("Selected rows contain duplicate SeriesInstanceUID values")
    if selected["trial_subject_id"].isna().any():
        raise ValueError("Selected series must have verified trial subject IDs")
    if not selected["visit_stage"].isin({"T0", "T1", "T2", "T3"}).all():
        raise ValueError("Selected series must have a verified T0/T1/T2/T3 stage")
    notes = selected["review_note"].astype("string").str.strip()
    if notes.isna().any() or notes.eq("").any() or notes.eq(DEFAULT_REVIEW_NOTE).any():
        raise ValueError("Replace the default review_note with a human review rationale")
    manifest = pd.DataFrame(
        {"SeriesInstanceUID": selected["series_instance_uid"].astype(str)}
    )
    summary = {
        "selected_series_count": int(len(manifest)),
        "selected_patient_count": int(selected["trial_subject_id"].nunique()),
        "selected_study_count": int(selected["study_instance_uid"].nunique())
        if "study_instance_uid" in selected.columns else None,
        "selected_visit_stage_counts": (
            selected["visit_stage"].value_counts().sort_index().astype(int).to_dict()
        ),
        "selected_modality_counts": (
            selected["modality"].value_counts().sort_index().astype(int).to_dict()
        ),
        "output_schema": ["SeriesInstanceUID"],
        "human_review_required": True,
        "privacy": "Manifest contains SeriesInstanceUIDs; keep local and private",
    }
    return manifest, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-csv", required=True)
    parser.add_argument("--manifest-csv", required=True)
    parser.add_argument("--summary-json", required=True)
    args = parser.parse_args()
    candidates = pd.read_csv(args.candidate_csv, dtype="string")
    manifest, summary = export_reviewed_manifest(candidates)
    output = Path(args.manifest_csv)
    summary_path = Path(args.summary_json)
    output.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(output, index=False)
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
