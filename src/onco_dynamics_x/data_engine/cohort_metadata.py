"""Normalise cBioPortal sample objects for safe cohort linkage."""
from __future__ import annotations

from typing import Any


def parse_cbioportal_sample_metadata(
    samples: list[dict[str, Any]],
) -> dict[str, dict[str, str]]:
    """Return sample ID -> patient ID and tissue type from cBioPortal Sample objects.

    cBioPortal Sample objects use id and type and may embed the Patient object
    in patient. Older client payloads can also expose sampleId and patientId
    directly; those aliases are accepted.
    """
    metadata: dict[str, dict[str, str]] = {}
    for sample in samples:
        sample_id = str(sample.get("sampleId") or sample.get("id") or "").strip()
        patient_object = sample.get("patient")
        patient = patient_object if isinstance(patient_object, dict) else {}
        patient_id = str(
            sample.get("patientId")
            or patient.get("patientId")
            or patient.get("id")
            or ""
        ).strip()
        sample_type = str(
            sample.get("type")
            or sample.get("sampleType")
            or sample.get("sample_type")
            or ""
        ).strip()
        if sample_id:
            metadata[sample_id] = {
                "patientId": patient_id,
                "type": sample_type,
            }
    return metadata
