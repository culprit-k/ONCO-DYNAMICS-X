"""Aggregate-only audit for biological calibration endpoint availability."""
from __future__ import annotations

import re
from typing import Any

from onco_dynamics_x.data_engine.cbioportal import CBioPortalClient, audit_cbioportal_study
from onco_dynamics_x.data_engine.molecular_audit import DEFAULT_COHORTS


def _normalise(value: Any) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(value or "").upper())


def _matches_target(target: str, attribute_id: Any, display_name: Any) -> bool:
    """Use specific clinical-attribute names; avoid substring false positives."""
    field_id = _normalise(attribute_id)
    label = _normalise(display_name)

    if target == "ki67_or_proliferation_index":
        return (
            field_id.startswith(("KI67", "MIB1", "PROLIFERATIONINDEX"))
            or label.startswith(("KI67", "MIB1", "PROLIFERATIONINDEX"))
        )
    if target == "er_status":
        return (
            field_id.startswith(("ERSTATUS", "ERIHC"))
            or label.startswith(("ERSTATUS", "ESTROGENRECEPTOR"))
        )
    if target == "pr_status":
        return (
            field_id.startswith(("PRSTATUS", "PRICH"))
            or label.startswith(("PRSTATUS", "PROGESTERONERECEPTOR"))
        )
    if target == "her2_status":
        return (
            field_id.startswith(("HER2STATUS", "HER2IHC", "HER2SNP6", "ERBB2STATUS"))
            or label.startswith(("HER2STATUS", "ERBB2STATUS"))
        )
    if target == "tumour_grade":
        return (
            field_id in {"GRADE", "TUMORGRADE", "TUMOURGRADE", "HISTOLOGICGRADE"}
            or label.endswith(("HISTOLOGICGRADE", "TUMORGRADE", "TUMOURGRADE"))
        )
    if target == "tumour_size":
        return (
            field_id.startswith(("TUMORSIZE", "TUMOURSIZE", "TUMORDIAMETER", "TUMOURDIAMETER"))
            or label.startswith(("TUMORSIZE", "TUMOURSIZE", "TUMORDIAMETER", "TUMOURDIAMETER"))
        )
    if target == "cellularity_or_hypoxia":
        target_prefixes = (
            "CELLULARITY",
            "BUFFAHYPOXIA",
            "RAGNUMHYPOXIA",
            "WINTERHYPOXIA",
            "HYPOXIA",
        )
        return field_id.startswith(target_prefixes) or label.startswith(target_prefixes)
    raise ValueError(f"Unknown calibration target category: {target}")


_TARGETS = (
    "ki67_or_proliferation_index",
    "er_status",
    "pr_status",
    "her2_status",
    "tumour_grade",
    "tumour_size",
    "cellularity_or_hypoxia",
)


def summarise_calibration_targets(
    study_id: str,
    study_name: str | None,
    attribute_coverage: list[dict[str, Any]],
) -> dict[str, Any]:
    """Match candidate target attributes by metadata name and report coverage."""
    targets: dict[str, list[dict[str, Any]]] = {}
    for target in _TARGETS:
        matches = []
        for row in attribute_coverage:
            if _matches_target(target, row.get("attribute_id"), row.get("display_name")):
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
