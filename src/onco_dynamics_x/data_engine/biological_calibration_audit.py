"""Aggregate-only audit for biological calibration endpoint availability."""
from __future__ import annotations

import re
from statistics import median
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
        prefixes = (
            "CELLULARITY",
            "BUFFAHYPOXIA",
            "RAGNUMHYPOXIA",
            "WINTERHYPOXIA",
            "HYPOXIA",
        )
        return field_id.startswith(prefixes) or label.startswith(prefixes)
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


def _nonempty(value: Any) -> bool:
    return value is not None and str(value).strip() not in {"", "NA", "NaN", "null"}


def summarise_target_value_overlap(
    attributes: list[dict[str, Any]],
    patient_records: list[dict[str, Any]],
    sample_records: list[dict[str, Any]],
    samples: list[dict[str, Any]],
    expression_sample_ids: set[str],
) -> list[dict[str, Any]]:
    """Summarise target value distributions and exact expression-list overlap.

    Entity identifiers are used only in memory for matching and are never
    returned. Ambiguous entities with conflicting values are counted separately.
    """
    sample_to_patient = {
        str(row["sampleId"]): str(row["patientId"])
        for row in samples
        if row.get("sampleId") and row.get("patientId")
    }
    expression_patients = {
        sample_to_patient[sample_id]
        for sample_id in expression_sample_ids
        if sample_id in sample_to_patient
    }
    output = []
    for target in _TARGETS:
        for attribute in attributes:
            attribute_id = str(attribute.get("clinicalAttributeId") or "")
            display = attribute.get("displayName") or attribute_id
            if not _matches_target(target, attribute_id, display):
                continue
            patient_level = bool(attribute.get("patientAttribute"))
            records = patient_records if patient_level else sample_records
            entity_key = "patientId" if patient_level else "sampleId"
            eligible = expression_patients if patient_level else expression_sample_ids
            grouped: dict[str, set[str]] = {}
            for row in records:
                if row.get("clinicalAttributeId") != attribute_id:
                    continue
                entity_id = row.get(entity_key)
                value = row.get("value")
                if entity_id and _nonempty(value):
                    grouped.setdefault(str(entity_id), set()).add(str(value).strip())
            relevant = set(grouped) & eligible
            unambiguous = {
                entity_id: next(iter(grouped[entity_id]))
                for entity_id in relevant
                if len(grouped[entity_id]) == 1
            }
            ambiguous_count = sum(
                len(grouped[entity_id]) > 1 for entity_id in relevant
            )
            counts: dict[str, int] = {}
            numeric_values: list[float] = []
            all_numeric = bool(unambiguous)
            for value in unambiguous.values():
                counts[value] = counts.get(value, 0) + 1
                try:
                    numeric_values.append(float(value))
                except ValueError:
                    all_numeric = False
            value_summary: dict[str, Any]
            if all_numeric and len(set(numeric_values)) > 20:
                ordered = sorted(numeric_values)
                value_summary = {
                    "numeric_min": min(ordered),
                    "numeric_median": float(median(ordered)),
                    "numeric_max": max(ordered),
                }
            else:
                value_summary = {
                    "category_counts": dict(sorted(counts.items())[:40]),
                    "categories_truncated": len(counts) > 40,
                }
            output.append(
                {
                    "target": target,
                    "attribute_id": attribute_id,
                    "display_name": display,
                    "level": "PATIENT" if patient_level else "SAMPLE",
                    "expression_entities_with_value": len(unambiguous),
                    "expression_entities_ambiguous": int(ambiguous_count),
                    "expression_entities_without_unambiguous_value": (
                        len(eligible) - len(unambiguous)
                    ),
                    "value_summary": value_summary,
                }
            )
    return output


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
            "method, units and temporal ordering must be verified before calibration."
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
        sample_list = client.get_json(f"sample-lists/{cohort.sample_list_id}")
        samples = client.get_all(
            f"studies/{cohort.study_id}/samples", {"projection": "SUMMARY"}
        )
        attributes = client.get_all(
            f"studies/{cohort.study_id}/clinical-attributes"
        )
        patient_records = client.get_all(
            f"studies/{cohort.study_id}/clinical-data",
            {"clinicalDataType": "PATIENT"},
        )
        sample_records = client.get_all(
            f"studies/{cohort.study_id}/clinical-data",
            {"clinicalDataType": "SAMPLE"},
        )
        summary = summarise_calibration_targets(
            cohort.study_id,
            report.get("study", {}).get("name"),
            report.get("clinical_inventory", {}).get("attribute_coverage", []),
        )
        expression_sample_ids = {
            str(value) for value in sample_list.get("sampleIds", [])
        }
        summary["expression_list_sample_count"] = len(expression_sample_ids)
        summary["target_value_overlap"] = summarise_target_value_overlap(
            attributes,
            patient_records,
            sample_records,
            samples,
            expression_sample_ids,
        )
        cohorts.append(summary)
    return {
        "analysis": "public_cBioPortal_breast_cancer_biological_calibration_target_audit",
        "cohorts": cohorts,
        "privacy": (
            "Aggregate metadata, category distributions and expression-list overlap only; "
            "no patient/sample identifiers, row-level clinical values, or expression "
            "values are exported."
        ),
        "calibration_readiness": {
            "longitudinal_tumour_volume_available_from_this_audit": False,
            "patient_matched_mri_cellularity_maps_available_from_this_audit": False,
            "serial_treatment_exposure_timeline_available_from_this_audit": False,
            "status": "cross_sectional_target_overlap_audit_only",
        },
    }
