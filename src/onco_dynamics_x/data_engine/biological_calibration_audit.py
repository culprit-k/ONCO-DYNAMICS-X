"""Aggregate-only audit for biological calibration endpoint availability."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

from onco_dynamics_x.data_engine.cbioportal import (
    CBioPortalClient,
    audit_cbioportal_study,
)
from onco_dynamics_x.data_engine.molecular_audit import DEFAULT_COHORTS


_TARGET_PATTERNS = {
    "ki67_or_proliferation_index": ("KI67", "MIB1", "PROLIFERATIONINDEX"),
    "er_status": ("ERSTATUS", "ESTROGENRECEPTOR"),
    "pr_status": ("PRSTATUS", "PROGESTERONERECEPTOR"),
    "her2_status": ("HER2STATUS", "ERBB2STATUS", "HER2IHC"),
    "tumour_grade": ("TUMORGRADE", "TUMOURGRADE", "HISTOLOGICGRADE", "GRADE"),
    "tumour_size": ("TUMORSIZE", "TUMOURSIZE", "TUMORDIAMETER", "TUMOURDIAMETER"),
    "cellularity_or_hypoxia": (
        "CELLULARITY", "ADC", "HYPOXIA", "OXYGEN", "HIF1",
    ),
}


def _normalise(value: Any) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(value or "").upper())


def summarise_calibration_targets(
    study_id: str,
    study_name: str | None,
    attribute_coverage: list[dict[str, Any]],
) -> dict[str, Any]:
    """Match candidate target attributes by metadata name and report coverage."""
    targets = {}
    for target, patterns in _TARGET_PATTERNS.items():
        matches = []
        for row in attribute_coverage:
            token = _normalise(
                f"{row.get('attribute_id', '')} {row.get('display_name', '')}"
            )
            if any(pattern in token for pattern in patterns):
                matches.append(
                    {
                        "attribute_id": row.get("attribute_id"),
                        "display_name": row.get("display_name"),
                        "level": row.get("level"),
                        "data_type": row.get("data_type"),
                        "entities_with_value": row.get("unique_entities_with_value"),
                        "records_with_value": row.get("records_with_value"),
                    }
                )
        targets[target] = matches
    return {
        "study_id": study_id,
        "study_name": study_name,
        "target_attribute_matches": targets,
        "matched_target_categories": [
            key for key, value in targets.items() if value
        ],
        "note": (
            "Attribute-name matches are candidates only. Source definitions, assay "
            "method, temporal ordering, units and sample correspondence must be "
            "verified before values are used as calibration targets."
        ),
    }


def run_audit(timeout: float = 60.0, retries: int = 3) -> dict[str, Any]:
    client = CBioPortalClient(timeout=timeout, retries=retries)
    cohorts = []
    for cohort in DEFAULT_COHORTS:
        report = audit_cbioportal_study(
            client=client,
            study_id=cohort.study_id,
            expression_sample_list_id=cohort.sample_list_id,
        )
        cohorts.append(
            summarise_calibration_targets(
                cohort.study_id,
                report.get("study", {}).get("name"),
                report.get("clinical_inventory", {}).get("attribute_coverage", []),
            )
        )
    return {
        "analysis": "public_cBioPortal_breast_cancer_biological_calibration_target_audit",
        "cohorts": cohorts,
        "privacy": (
            "Aggregate attribute metadata and entity counts only; no patient/sample "
            "identifiers, row-level clinical values, or expression values are exported."
        ),
        "calibration_readiness": {
            "longitudinal_tumour_volume_available_from_this_audit": False,
            "patient_matched_mri_cellularity_maps_available_from_this_audit": False,
            "serial_treatment_exposure_timeline_available_from_this_audit": False,
            "status": "metadata_feasibility_audit_only",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Audit public TCGA/METABRIC clinical measurement availability."
    )
    parser.add_argument(
        "--output",
        default="results/biological_calibration/target_availability.json",
    )
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--retries", type=int, default=3)
    args = parser.parse_args()
    report = run_audit(timeout=args.timeout, retries=args.retries)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for cohort in report["cohorts"]:
        print(
            f"{cohort['study_id']}: candidate target categories "
            f"{cohort['matched_target_categories']}"
        )
    print(f"Aggregate audit: {output}")


if __name__ == "__main__":
    main()
