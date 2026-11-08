"""Run an exploratory age+stage Cox model trained in TCGA and tested in METABRIC."""
from __future__ import annotations

import argparse

from onco_dynamics_x.data_engine.cbioportal import CBioPortalClient
from onco_dynamics_x.data_engine.survival_analysis import (
    fetch_survival_cohort,
    run_cross_cohort_cox,
    write_survival_report,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results/clinical_survival_baseline/report.json")
    parser.add_argument("--penalizer", type=float, default=0.1)
    args = parser.parse_args()

    client = CBioPortalClient(timeout=60, retries=3)
    tcga, tcga_qc = fetch_survival_cohort(
        client=client,
        study_id="brca_tcga_pan_can_atlas_2018",
        sample_list_id="brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
        age_attribute="AGE",
        stage_attribute="AJCC_PATHOLOGIC_TUMOR_STAGE",
        stage_level="PATIENT",
    )
    metabric, metabric_qc = fetch_survival_cohort(
        client=client,
        study_id="brca_metabric",
        sample_list_id="brca_metabric_mrna",
        age_attribute="AGE_AT_DIAGNOSIS",
        stage_attribute="TUMOR_STAGE",
        stage_level="SAMPLE",
    )
    model_report = run_cross_cohort_cox(tcga, metabric, penalizer=args.penalizer)
    report = {
        "tcga_data_quality": tcga_qc,
        "metabric_data_quality": metabric_qc,
        "model": model_report,
        "data_governance": {
            "raw_patient_rows_uploaded": False,
            "sample_ids_uploaded": False,
            "only_aggregate_metrics_uploaded": True,
        },
    }
    write_survival_report(report, args.output)
    print("TCGA patients retained:", tcga_qc["patients_retained_for_model"])
    print("METABRIC patients retained:", metabric_qc["patients_retained_for_model"])
    print("External METABRIC Harrell C-index:", model_report["external_cohort"]["harrell_c_index"])
    print("Aggregate report:", args.output)


if __name__ == "__main__":
    main()
