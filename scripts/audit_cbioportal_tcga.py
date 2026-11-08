"""Audit TCGA-BRCA PanCancer Atlas clinical and molecular metadata."""
from __future__ import annotations

import argparse

from onco_dynamics_x.data_engine.cbioportal import (
    CBioPortalClient,
    audit_cbioportal_study,
    write_audit_outputs,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Audit the TCGA-BRCA PanCancer Atlas development cohort."
    )
    parser.add_argument("--output-json", default="results/tcga_cbioportal_audit/report.json")
    parser.add_argument("--output-csv", default="results/tcga_cbioportal_audit/clinical_attribute_coverage.csv")
    args = parser.parse_args()
    report = audit_cbioportal_study(
        client=CBioPortalClient(),
        study_id="brca_tcga_pan_can_atlas_2018",
    )
    write_audit_outputs(report, args.output_json, args.output_csv)
    print(f"Study: {report['study'].get('name')}")
    print(f"Patients: {report['cohort_counts']['unique_patients']}")
    print(f"Samples: {report['cohort_counts']['unique_samples']}")
    print(f"Molecular profiles: {len(report['molecular_profiles'])}")
    print(f"Clinical attributes: {report['clinical_inventory']['clinical_attribute_count']}")
    print(f"Aggregate report: {args.output_json}")


if __name__ == "__main__":
    main()
