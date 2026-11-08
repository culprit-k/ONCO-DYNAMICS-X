"""Read-only cBioPortal public API client and aggregate cohort audit.

The audit deliberately writes aggregate metadata only. It never exports
patient/sample identifiers or individual clinical values.
"""
from __future__ import annotations

import json
import time
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


class CBioPortalError(RuntimeError):
    """Raised when the cBioPortal API cannot provide a valid response."""


class CBioPortalClient:
    """Small standard-library client for public cBioPortal API v3 endpoints."""

    def __init__(self, base_url="https://www.cbioportal.org/api/v3", timeout=45.0,
                 retries=3, retry_delay=1.0):
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.retries = retries
        self.retry_delay = retry_delay

    def get_json(self, path: str, params: dict[str, Any] | None = None) -> Any:
        query = f"?{urlencode(params)}" if params else ""
        url = f"{self.base_url}/{path.lstrip('/')}{query}"
        request = Request(url, headers={
            "Accept": "application/json",
            "User-Agent": "ONCO-DYNAMICS-X-cohort-audit/1.0",
        })
        last_error = None
        for attempt in range(self.retries + 1):
            try:
                with urlopen(request, timeout=self.timeout) as response:
                    payload = response.read().decode("utf-8")
                return json.loads(payload)
            except HTTPError as exc:
                last_error = exc
                if exc.code not in {429, 500, 502, 503, 504} or attempt >= self.retries:
                    detail = exc.read().decode("utf-8", errors="replace")[:500]
                    raise CBioPortalError(
                        f"cBioPortal returned HTTP {exc.code} for {path}: {detail}"
                    ) from exc
            except (URLError, TimeoutError, json.JSONDecodeError) as exc:
                last_error = exc
                if attempt >= self.retries:
                    break
            time.sleep(self.retry_delay * (2**attempt))
        raise CBioPortalError(f"Could not retrieve {path}: {last_error}") from last_error

    def get_all(self, path: str, params: dict[str, Any] | None = None,
                page_size: int = 1000) -> list[dict[str, Any]]:
        """Retrieve a paginated list endpoint, stopping on a short page."""
        base_params = dict(params or {})
        records = []
        page_number = 0
        while True:
            page_params = {**base_params, "pageSize": page_size, "pageNumber": page_number}
            page = self.get_json(path, page_params)
            if not isinstance(page, list):
                raise CBioPortalError(
                    f"Expected a JSON list from {path}; received {type(page).__name__}."
                )
            records.extend(item for item in page if isinstance(item, dict))
            if len(page) < page_size:
                return records
            page_number += 1


def _nonempty(value: Any) -> bool:
    return value is not None and str(value).strip() not in {"", "NA", "NaN", "null"}


def _normalise(text: Any) -> str:
    return "".join(char for char in str(text or "").upper() if char.isalnum())


def _candidate_endpoint(attribute: dict[str, Any]) -> str | None:
    """Classify endpoint fields by name only; classification is not validation."""
    token = _normalise(
        f"{attribute.get('clinicalAttributeId', '')} "
        f"{attribute.get('displayName', '')} {attribute.get('description', '')}"
    )
    if any(key in token for key in ("OVERALLSURVIVAL", "OSSTATUS", "OSMONTHS", "OSYEARS")):
        return "overall_survival_candidate"
    if any(key in token for key in ("DISEASESPECIFICSURVIVAL", "DSSSTATUS", "DSSMONTHS", "DSSYEARS")):
        return "disease_specific_survival_candidate"
    if any(key in token for key in ("RELAPSEFREESURVIVAL", "RELAPSEFREE", "RFSSTATUS", "RFSMONTHS")):
        return "relapse_free_candidate"
    if any(key in token for key in ("DISEASEFREE", "DFSSTATUS", "DFSMONTHS")):
        return "disease_free_candidate"
    if any(key in token for key in ("PROGRESSIONFREESURVIVAL", "PFSSTATUS", "PFSMONTHS")):
        return "progression_free_candidate"
    return None


def _candidate_treatment_field(attribute: dict[str, Any]) -> bool:
    token = _normalise(
        f"{attribute.get('clinicalAttributeId', '')} "
        f"{attribute.get('displayName', '')} {attribute.get('description', '')}"
    )
    return any(term in token for term in (
        "TREATMENT", "CHEMOTHERAPY", "RADIOTHERAPY", "RADIATION",
        "ENDOCRINE", "HORMONETHERAPY", "TARGETEDTHERAPY",
        "TREATMENTTYPE", "THERAPYTYPE",
    ))


def summarize_clinical_coverage(attributes, patient_data, sample_data):
    """Return per-attribute aggregate coverage without exposing entity IDs."""
    summaries = []
    for attribute in attributes:
        attribute_id = attribute.get("clinicalAttributeId")
        if not attribute_id:
            continue
        is_patient = bool(attribute.get("patientAttribute"))
        rows = patient_data if is_patient else sample_data
        matching = [row for row in rows if row.get("clinicalAttributeId") == attribute_id]
        populated = [row for row in matching if _nonempty(row.get("value"))]
        entity_key = "patientId" if is_patient else "sampleId"
        summaries.append({
            "attribute_id": attribute_id,
            "display_name": attribute.get("displayName") or attribute.get("description") or attribute_id,
            "data_type": attribute.get("datatype") or attribute.get("dataType"),
            "level": "PATIENT" if is_patient else "SAMPLE",
            "records_with_value": len(populated),
            "unique_entities_with_value": len({
                row.get(entity_key) for row in populated if row.get(entity_key)
            }),
            "candidate_survival_endpoint": _candidate_endpoint(attribute),
            "candidate_treatment_field": _candidate_treatment_field(attribute),
        })
    return sorted(summaries, key=lambda item: (item["level"], item["attribute_id"]))


def audit_metabric(client: CBioPortalClient | None = None) -> dict[str, Any]:
    """Audit METABRIC study metadata, molecular profiles and clinical coverage."""
    api = client or CBioPortalClient()
    study_id = "brca_metabric"
    study = api.get_json(f"studies/{study_id}")
    patients = api.get_all(f"studies/{study_id}/patients", {"projection": "SUMMARY"})
    samples = api.get_all(f"studies/{study_id}/samples", {"projection": "SUMMARY"})
    attributes = api.get_all(f"studies/{study_id}/clinical-attributes")
    profiles = api.get_all(f"studies/{study_id}/molecular-profiles")
    sample_lists = api.get_all(f"studies/{study_id}/sample-lists")
    patient_clinical = api.get_all(
        f"studies/{study_id}/clinical-data", {"clinicalDataType": "PATIENT"}
    )
    sample_clinical = api.get_all(
        f"studies/{study_id}/clinical-data", {"clinicalDataType": "SAMPLE"}
    )

    patient_ids = {item.get("patientId") for item in patients if item.get("patientId")}
    sample_ids = {item.get("sampleId") for item in samples if item.get("sampleId")}
    sample_patient_ids = {item.get("patientId") for item in samples if item.get("patientId")}
    sample_patient_links = sum(
        1 for item in samples
        if item.get("sampleId") and item.get("patientId") in patient_ids
    )
    orphan_samples = len(samples) - sample_patient_links
    coverage = summarize_clinical_coverage(attributes, patient_clinical, sample_clinical)
    survival_candidates = [row for row in coverage if row["candidate_survival_endpoint"]]
    treatment_candidates = [row for row in coverage if row["candidate_treatment_field"]]
    profile_summaries = [{
        "molecular_profile_id": profile.get("molecularProfileId"),
        "name": profile.get("name"),
        "description": profile.get("description"),
        "molecular_alteration_type": profile.get("molecularAlterationType"),
        "datatype": profile.get("datatype"),
        "patient_level": profile.get("patientLevel"),
        "generic_assay_type": profile.get("genericAssayType"),
    } for profile in profiles]
    patient_endpoint_candidates = [
        row for row in survival_candidates
        if row["level"] == "PATIENT" and row["unique_entities_with_value"] > 0
    ]
    patient_treatment_candidates = [
        row for row in treatment_candidates
        if row["level"] == "PATIENT" and row["unique_entities_with_value"] > 0
    ]

    return {
        "study": {
            "study_id": study_id,
            "name": study.get("name"),
            "description": study.get("description"),
            "cancer_type_id": study.get("cancerTypeId"),
            "type_of_cancer_id": study.get("typeOfCancerId"),
            "public_study": study.get("publicStudy"),
            "reference_genome": study.get("referenceGenome"),
            "citation": study.get("citation"),
            "pmid": study.get("pmid"),
        },
        "cohort_counts": {
            "unique_patients": len(patient_ids),
            "unique_samples": len(sample_ids),
            "patients_with_at_least_one_sample": len(sample_patient_ids & patient_ids),
            "samples_linked_to_known_patient": sample_patient_links,
            "samples_without_patient_link": orphan_samples,
        },
        "molecular_profiles": profile_summaries,
        "sample_lists": [{
            "sample_list_id": item.get("sampleListId"),
            "name": item.get("name"),
            "description": item.get("description"),
            "category": item.get("category"),
        } for item in sample_lists],
        "clinical_inventory": {
            "clinical_attribute_count": len(attributes),
            "patient_clinical_records": len(patient_clinical),
            "sample_clinical_records": len(sample_clinical),
            "attribute_coverage": coverage,
            "patient_level_survival_field_candidates": patient_endpoint_candidates,
            "patient_level_treatment_field_candidates": patient_treatment_candidates,
        },
        "audit_checks": {
            "study_metadata_present": study.get("studyId") == study_id,
            "patients_have_unique_ids": len(patient_ids) == len(patients),
            "samples_have_unique_ids": len(sample_ids) == len(samples),
            "all_samples_link_to_known_patient": orphan_samples == 0,
            "has_molecular_profiles": bool(profiles),
            "has_patient_clinical_attributes": any(
                bool(item.get("patientAttribute")) for item in attributes
            ),
            "has_patient_survival_field_candidate": bool(patient_endpoint_candidates),
            "has_patient_treatment_field_candidate": bool(patient_treatment_candidates),
        },
        "interpretation": {
            "endpoint_candidates_are_name_based_only": True,
            "clinical_endpoint_definitions_verified": False,
            "survival_model_fit": False,
            "external_validation_performed": False,
            "treatment_response_validated": False,
            "raw_patient_or_sample_records_exported": False,
            "note": (
                "Candidate fields require data-dictionary and source-publication "
                "verification. Field presence alone does not establish endpoint "
                "definition, time origin, censoring, treatment exposure or response."
            ),
        },
    }


def write_audit_outputs(report: dict[str, Any], json_path: str, csv_path: str) -> None:
    """Write aggregate-only JSON and clinical-attribute coverage CSV."""
    import csv
    from pathlib import Path

    json_output = Path(json_path)
    csv_output = Path(csv_path)
    json_output.parent.mkdir(parents=True, exist_ok=True)
    csv_output.parent.mkdir(parents=True, exist_ok=True)
    json_output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    rows = report["clinical_inventory"]["attribute_coverage"]
    fields = [
        "attribute_id", "display_name", "data_type", "level", "records_with_value",
        "unique_entities_with_value", "candidate_survival_endpoint", "candidate_treatment_field",
    ]
    with csv_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
