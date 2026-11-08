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
    """Small standard-library client for public cBioPortal REST endpoints."""

    def __init__(
        self,
        base_url: str = "https://www.cbioportal.org/api",
        timeout: float = 45.0,
        retries: int = 3,
        retry_delay: float = 1.0,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.retries = retries
        self.retry_delay = retry_delay

    def get_json(self, path: str, params: dict[str, Any] | None = None) -> Any:
        query = f"?{urlencode(params)}" if params else ""
        url = f"{self.base_url}/{path.lstrip('/')}{query}"
        request = Request(
            url,
            headers={
                "Accept": "application/json",
                "User-Agent": "ONCO-DYNAMICS-X-cohort-audit/1.0",
            },
        )
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

    def get_all(
        self,
        path: str,
        params: dict[str, Any] | None = None,
        page_size: int = 1000,
    ) -> list[dict[str, Any]]:
        """Retrieve a paginated list endpoint, stopping on a short page."""
        base_params = dict(params or {})
        records = []
        page_number = 0
        while True:
            page_params = {
                **base_params,
                "pageSize": page_size,
                "pageNumber": page_number,
            }
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
    if any(
        key in token
        for key in ("DISEASESPECIFICSURVIVAL", "DSSSTATUS", "DSSMONTHS", "DSSYEARS")
    ):
        return "disease_specific_survival_candidate"
    if any(
        key in token
        for key in ("RELAPSEFREESURVIVAL", "RELAPSEFREE", "RFSSTATUS", "RFSMONTHS")
    ):
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
    return any(
        term in token
        for term in (
            "TREATMENT",
            "CHEMOTHERAPY",
            "RADIOTHERAPY",
            "RADIATION",
            "ENDOCRINE",
            "HORMONETHERAPY",
            "TARGETEDTHERAPY",
            "TREATMENTTYPE",
            "THERAPYTYPE",
        )
    )


def summarize_clinical_coverage(
    attributes: list[dict[str, Any]],
    patient_data: list[dict[str, Any]],
    sample_data: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Return per-attribute aggregate coverage without exposing entity IDs."""
    summaries = []
    for attribute in attributes:
        attribute_id = attribute.get("clinicalAttributeId")
        if not attribute_id:
            continue
        is_patient = bool(attribute.get("patientAttribute"))
        rows = patient_data if is_patient else sample_data
        matching = [
            row for row in rows if row.get("clinicalAttributeId") == attribute_id
        ]
        populated = [row for row in matching if _nonempty(row.get("value"))]
        entity_key = "patientId" if is_patient else "sampleId"
        display_name = (
            attribute.get("displayName")
            or attribute.get("description")
            or attribute_id
        )
        summaries.append(
            {
                "attribute_id": attribute_id,
                "display_name": display_name,
                "data_type": attribute.get("datatype") or attribute.get("dataType"),
                "level": "PATIENT" if is_patient else "SAMPLE",
                "records_with_value": len(populated),
                "unique_entities_with_value": len(
                    {
                        row.get(entity_key)
                        for row in populated
                        if row.get(entity_key)
                    }
                ),
                "candidate_survival_endpoint": _candidate_endpoint(attribute),
                "candidate_treatment_field": _candidate_treatment_field(attribute),
            }
        )
    return sorted(summaries, key=lambda item: (item["level"], item["attribute_id"]))


def _patient_attribute_values(records, attribute_id):
    grouped = {}
    for row in records:
        patient_id = row.get("patientId")
        value = row.get("value")
        if patient_id and row.get("clinicalAttributeId") == attribute_id and _nonempty(value):
            grouped.setdefault(patient_id, set()).add(str(value).strip())
    values = {}
    ambiguous = 0
    for patient_id, observed in grouped.items():
        if len(observed) == 1:
            values[patient_id] = next(iter(observed))
        else:
            ambiguous += 1
    return values, ambiguous


def summarize_survival_endpoint_pairs(patient_data):
    """Measure paired survival-field completeness without exporting patient IDs."""
    definitions = (
        ("overall_survival", "OS_MONTHS", "OS_STATUS"),
        ("relapse_free_survival", "RFS_MONTHS", "RFS_STATUS"),
    )
    summaries = []
    for endpoint, time_field, status_field in definitions:
        times, ambiguous_times = _patient_attribute_values(patient_data, time_field)
        statuses, ambiguous_statuses = _patient_attribute_values(
            patient_data, status_field
        )
        complete_ids = set(times) & set(statuses)
        parsed_times = []
        unparseable = 0
        for patient_id in complete_ids:
            try:
                parsed_times.append(float(times[patient_id]))
            except (TypeError, ValueError):
                unparseable += 1
        positive = [value for value in parsed_times if value > 0]
        status_counts = {}
        for value in statuses.values():
            status_counts[value] = status_counts.get(value, 0) + 1
        summaries.append(
            {
                "endpoint": endpoint,
                "time_field": time_field,
                "status_field": status_field,
                "patients_with_time": len(times),
                "patients_with_status": len(statuses),
                "patients_with_both_fields": len(complete_ids),
                "ambiguous_time_patient_count": ambiguous_times,
                "ambiguous_status_patient_count": ambiguous_statuses,
                "complete_pairs_with_numeric_time": len(parsed_times),
                "positive_time_count": len(positive),
                "zero_or_negative_time_count": len(
                    [value for value in parsed_times if value <= 0]
                ),
                "unparseable_time_count": unparseable,
                "time_range_months": (
                    [round(min(positive), 4), round(max(positive), 4)]
                    if positive
                    else None
                ),
                "status_value_counts": dict(sorted(status_counts.items())),
            }
        )
    return summaries


def summarize_expression_survival_linkage(
    samples, expression_sample_ids, patient_clinical
):
    """Summarize expression-sample overlap with complete patient survival fields."""
    sample_to_patient = {
        row.get("sampleId"): row.get("patientId")
        for row in samples
        if row.get("sampleId") and row.get("patientId")
    }
    expression_ids = set(expression_sample_ids)
    expression_patient_ids = {
        sample_to_patient[sample_id]
        for sample_id in expression_ids
        if sample_id in sample_to_patient
    }
    os_time, _ = _patient_attribute_values(patient_clinical, "OS_MONTHS")
    os_status, _ = _patient_attribute_values(patient_clinical, "OS_STATUS")
    rfs_time, _ = _patient_attribute_values(patient_clinical, "RFS_MONTHS")
    rfs_status, _ = _patient_attribute_values(patient_clinical, "RFS_STATUS")
    os_complete = set(os_time) & set(os_status)
    rfs_complete = set(rfs_time) & set(rfs_status)
    return {
        "expression_sample_list_sample_count": len(expression_ids),
        "expression_samples_linked_to_patient": sum(
            sample_id in sample_to_patient for sample_id in expression_ids
        ),
        "unique_expression_patients": len(expression_patient_ids),
        "expression_patients_with_complete_overall_survival": len(
            expression_patient_ids & os_complete
        ),
        "expression_patients_with_complete_relapse_free_survival": len(
            expression_patient_ids & rfs_complete
        ),
    }


def audit_cbioportal_study(
    client: CBioPortalClient | None = None,
    study_id: str = "brca_metabric",
    expression_sample_list_id: str | None = None,
) -> dict[str, Any]:
    """Audit a cBioPortal study without exporting row-level patient data."""
    api = client or CBioPortalClient()
    study = api.get_json(f"studies/{study_id}")
    patients = api.get_all(
        f"studies/{study_id}/patients", {"projection": "SUMMARY"}
    )
    samples = api.get_all(
        f"studies/{study_id}/samples", {"projection": "SUMMARY"}
    )
    attributes = api.get_all(f"studies/{study_id}/clinical-attributes")
    profiles = api.get_all(f"studies/{study_id}/molecular-profiles")
    sample_lists = api.get_all(f"studies/{study_id}/sample-lists")
    expression_sample_list = (
        api.get_json(f"sample-lists/{expression_sample_list_id}")
        if expression_sample_list_id
        else None
    )
    patient_clinical = api.get_all(
        f"studies/{study_id}/clinical-data", {"clinicalDataType": "PATIENT"}
    )
    sample_clinical = api.get_all(
        f"studies/{study_id}/clinical-data", {"clinicalDataType": "SAMPLE"}
    )

    patient_ids = {
        item.get("patientId") for item in patients if item.get("patientId")
    }
    sample_ids = {
        item.get("sampleId") for item in samples if item.get("sampleId")
    }
    sample_patient_ids = {
        item.get("patientId") for item in samples if item.get("patientId")
    }
    sample_patient_links = sum(
        1
        for item in samples
        if item.get("sampleId") and item.get("patientId") in patient_ids
    )
    orphan_samples = len(samples) - sample_patient_links
    coverage = summarize_clinical_coverage(
        attributes, patient_clinical, sample_clinical
    )
    survival_candidates = [
        row for row in coverage if row["candidate_survival_endpoint"]
    ]
    treatment_candidates = [
        row for row in coverage if row["candidate_treatment_field"]
    ]
    profile_summaries = [
        {
            "molecular_profile_id": profile.get("molecularProfileId"),
            "name": profile.get("name"),
            "description": profile.get("description"),
            "molecular_alteration_type": profile.get("molecularAlterationType"),
            "datatype": profile.get("datatype"),
            "patient_level": profile.get("patientLevel"),
            "generic_assay_type": profile.get("genericAssayType"),
        }
        for profile in profiles
    ]
    patient_endpoint_candidates = [
        row
        for row in survival_candidates
        if row["level"] == "PATIENT" and row["unique_entities_with_value"] > 0
    ]
    patient_treatment_candidates = [
        row
        for row in treatment_candidates
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
            "patients_with_at_least_one_sample": len(
                sample_patient_ids & patient_ids
            ),
            "samples_linked_to_known_patient": sample_patient_links,
            "samples_without_patient_link": orphan_samples,
        },
        "molecular_profiles": profile_summaries,
        "expression_survival_linkage": (
            summarize_expression_survival_linkage(
                samples, expression_sample_list.get("sampleIds", []), patient_clinical
            )
            if expression_sample_list is not None
            else None
        ),
        "sample_lists": [
            {
                "sample_list_id": item.get("sampleListId"),
                "name": item.get("name"),
                "description": item.get("description"),
                "category": item.get("category"),
            }
            for item in sample_lists
        ],
        "clinical_inventory": {
            "clinical_attribute_count": len(attributes),
            "patient_clinical_records": len(patient_clinical),
            "sample_clinical_records": len(sample_clinical),
            "attribute_coverage": coverage,
            "patient_level_survival_field_candidates": patient_endpoint_candidates,
            "patient_level_treatment_field_candidates": patient_treatment_candidates,
            "survival_endpoint_pair_summary": summarize_survival_endpoint_pairs(
                patient_clinical
            ),
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
            "has_patient_survival_field_candidate": bool(
                patient_endpoint_candidates
            ),
            "has_patient_treatment_field_candidate": bool(
                patient_treatment_candidates
            ),
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




def audit_metabric(client: CBioPortalClient | None = None) -> dict[str, Any]:
    """Audit the METABRIC study and its mRNA-expression sample-list overlap."""
    return audit_cbioportal_study(
        client=client,
        study_id="brca_metabric",
        expression_sample_list_id="brca_metabric_mrna",
    )


def write_audit_outputs(
    report: dict[str, Any], json_path: str, csv_path: str
) -> None:
    """Write aggregate-only JSON and clinical-attribute coverage CSV."""
    import csv
    from pathlib import Path

    json_output = Path(json_path)
    csv_output = Path(csv_path)
    json_output.parent.mkdir(parents=True, exist_ok=True)
    csv_output.parent.mkdir(parents=True, exist_ok=True)
    json_output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    rows = report["clinical_inventory"]["attribute_coverage"]
    fields = [
        "attribute_id",
        "display_name",
        "data_type",
        "level",
        "records_with_value",
        "unique_entities_with_value",
        "candidate_survival_endpoint",
        "candidate_treatment_field",
    ]
    with csv_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
