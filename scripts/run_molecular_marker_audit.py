"""Run the aggregate-only TCGA-BRCA / METABRIC marker availability audit."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from onco_dynamics_x.data_engine.molecular_audit import (
    audit_breast_cancer_marker_panel,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default="results/molecular_audit/marker_panel_overlap.json",
    )
    args = parser.parse_args()
    report = audit_breast_cancer_marker_panel()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("Molecular marker availability audit")
    for cohort in report["cohorts"]:
        print(
            cohort["cohort_label"],
            "samples:", cohort["sample_list_sample_count"],
            "markers with finite values:", cohort["markers_with_finite_values"],
            "/", cohort["marker_panel_size"],
        )
    print(
        "Markers present in all queried profiles:",
        report["cross_cohort_overlap"]["marker_count_present_in_all_profiles"],
    )
    print("Aggregate-only report:", output)


if __name__ == "__main__":
    main()
