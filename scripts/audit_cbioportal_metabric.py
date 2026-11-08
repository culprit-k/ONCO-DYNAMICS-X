"""Run a read-only aggregate audit of METABRIC through cBioPortal."""
from __future__ import annotations

import argparse

from onco_dynamics_x.data_engine.cbioportal import (
    CBioPortalClient,
    audit_metabric,
    write_audit_outputs,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Audit public METABRIC cohort metadata and clinical field coverage."
    )
    parser.add_argument("--output-json", default="results/metabric_audit/report.json")
    parser.add_argument("--output-csv", default="results/metabric_audit/clinical_attribute_coverage.csv")
    parser.add_argument("--timeout", type=float, default=45.0)
    parser.add_argument("--retries", type=int, default=3)
    args = parser.parse_args()

    report = audit_metabric(CBioPortalClient(timeout=args.timeout, retries=args.retries))
    write_audit_outputs(report, args.output_json, args.output_csv)

    print(f"Study: {report['study'].get('name')}")
    print(f"Patients: {report['cohort_counts']['unique_patients']}")
    print(f"Samples: {report['cohort_counts']['unique_samples']}")
    print(f"Molecular profiles: {len(report['molecular_profiles'])}")
    print(f"Clinical attributes: {report['clinical_inventory']['clinical_attribute_count']}")
    print("Patient-level survival field candidates: "
          f"{len(report['clinical_inventory']['patient_level_survival_field_candidates'])}")
    print("Patient-level treatment field candidates: "
          f"{len(report['clinical_inventory']['patient_level_treatment_field_candidates'])}")
    print(f"Aggregate report: {args.output_json}")
    print(f"Attribute inventory: {args.output_csv}")


if __name__ == "__main__":
    main()
