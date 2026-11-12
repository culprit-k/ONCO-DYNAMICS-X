"""Create a minimal TCIA Data Retriever CSV from a reviewed series selection.

The patient/visit-to-series choice is explicit and human-reviewed. The output
contains SeriesInstanceUID only; do not commit selection/inventory files.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

REQUIRED_INVENTORY = {"patient_id", "visit_stage", "series_instance_uid"}
REQUIRED_SELECTION = {
    "patient_id", "visit_stage", "series_instance_uid", "include",
}


def build_series_manifest(
    inventory: pd.DataFrame,
    selection: pd.DataFrame,
) -> tuple[pd.DataFrame, dict]:
    missing_inventory = REQUIRED_INVENTORY - set(inventory.columns)
    missing_selection = REQUIRED_SELECTION - set(selection.columns)
    if missing_inventory or missing_selection:
        raise ValueError(
            f"Missing columns; inventory={sorted(missing_inventory)}, "
            f"selection={sorted(missing_selection)}"
        )
    inv = inventory.copy()
    sel = selection.copy()
    for frame in (inv, sel):
        frame["patient_id"] = frame["patient_id"].astype("string").str.strip().str.zfill(6)
        frame["visit_stage"] = frame["visit_stage"].astype("string").str.strip().str.upper()
        frame["series_instance_uid"] = frame["series_instance_uid"].astype("string").str.strip()
    if inv["series_instance_uid"].duplicated().any():
        raise ValueError("Inventory contains duplicate SeriesInstanceUID values")
    if sel["series_instance_uid"].duplicated().any():
        raise ValueError("Selection contains duplicate SeriesInstanceUID values")
    selected_flags = sel["include"].astype(str).str.strip().str.lower()
    valid_flags = {"1", "0", "true", "false", "yes", "no"}
    if not set(selected_flags.unique()).issubset(valid_flags):
        raise ValueError("include must use 1/0, true/false or yes/no")
    include = selected_flags.isin({"1", "true", "yes"})
    chosen = sel.loc[include].copy()
    if chosen.empty:
        raise ValueError("No series marked include=1")

    if chosen[["patient_id", "visit_stage", "series_instance_uid"]].isna().any().any():
        raise ValueError("Selected series must have patient ID, visit stage and series UID")
    if not chosen["visit_stage"].isin({"T0", "T1", "T2", "T3"}).all():
        raise ValueError("visit_stage must be one of T0, T1, T2, T3")

    check = chosen.merge(
        inv[["patient_id", "visit_stage", "series_instance_uid"]],
        on=["patient_id", "visit_stage", "series_instance_uid"],
        how="left",
        validate="one_to_one",
        indicator=True,
    )
    if not check["_merge"].eq("both").all():
        invalid = int((check["_merge"] != "both").sum())
        raise ValueError(f"{invalid} selected series do not match the local inventory")

    manifest = pd.DataFrame({"SeriesInstanceUID": chosen["series_instance_uid"].astype(str)})
    manifest = manifest.drop_duplicates().reset_index(drop=True)
    summary = {
        "selected_series_count": int(len(manifest)),
        "selected_patient_count": int(chosen["patient_id"].nunique()),
        "selected_visit_counts": chosen.groupby("visit_stage").size().astype(int).to_dict(),
        "output_schema": ["SeriesInstanceUID"],
        "download_policy": "only explicitly selected series; no full-collection download",
        "privacy": "Selection and inventory files contain patient IDs; manifest contains series UIDs and must remain local.",
    }
    return manifest, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory-csv", required=True, help="Private local DICOM series inventory")
    parser.add_argument("--selection-csv", required=True, help="Private reviewed patient/visit/series selection")
    parser.add_argument("--manifest-csv", required=True, help="Local SeriesInstanceUID manifest")
    parser.add_argument("--summary-json", required=True, help="Aggregate-only summary")
    args = parser.parse_args()
    inventory = pd.read_csv(args.inventory_csv, dtype="string")
    selection = pd.read_csv(args.selection_csv, dtype="string")
    manifest, summary = build_series_manifest(inventory, selection)
    output = Path(args.manifest_csv)
    summary_path = Path(args.summary_json)
    output.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(output, index=False)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
