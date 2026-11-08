"""Cohort manifest validation without reading raw participant data."""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

_REQUIRED_COLUMNS = {
    "cohort_id",
    "patient_id",
    "specimen_id",
    "sample_id",
    "assay_type",
    "modality",
    "data_access",
    "source_accession",
}
_ALLOWED_ACCESS = {"open", "controlled", "restricted", "unknown"}


@dataclass(frozen=True)
class ManifestAudit:
    """Machine-readable summary of a sample metadata manifest."""

    rows: int
    unique_patients: int
    unique_samples: int
    errors: tuple[str, ...]
    warnings: tuple[str, ...]

    @property
    def passed(self) -> bool:
        return not self.errors


def audit_manifest(path: str | Path) -> ManifestAudit:
    """Validate required columns, identifiers and data-access labels in a CSV."""
    manifest_path = Path(path)
    if not manifest_path.exists():
        raise FileNotFoundError(manifest_path)
    with manifest_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        columns = set(reader.fieldnames or [])
        missing = sorted(_REQUIRED_COLUMNS - columns)
        if missing:
            return ManifestAudit(
                rows=0,
                unique_patients=0,
                unique_samples=0,
                errors=(f"Missing required columns: {', '.join(missing)}",),
                warnings=(),
            )
        records = list(reader)

    errors: list[str] = []
    warnings: list[str] = []
    samples: set[str] = set()
    patients: set[str] = set()
    patient_by_sample: dict[str, str] = {}

    for line_number, row in enumerate(records, start=2):
        for field in _REQUIRED_COLUMNS:
            if not (row.get(field) or "").strip():
                errors.append(f"Line {line_number}: '{field}' is empty.")
        sample_id = (row.get("sample_id") or "").strip()
        patient_id = (row.get("patient_id") or "").strip()
        access = (row.get("data_access") or "").strip().lower()
        if sample_id:
            if sample_id in samples:
                errors.append(f"Line {line_number}: duplicate sample_id '{sample_id}'.")
            samples.add(sample_id)
            patient_by_sample[sample_id] = patient_id
        if patient_id:
            patients.add(patient_id)
        if access and access not in _ALLOWED_ACCESS:
            errors.append(
                f"Line {line_number}: unsupported data_access '{access}'."
            )
        if access in {"controlled", "restricted", "unknown"}:
            warnings.append(
                f"Line {line_number}: {access} data must not be uploaded to public Git."
            )

    if not records:
        warnings.append("Manifest contains no sample rows.")
    return ManifestAudit(
        rows=len(records),
        unique_patients=len(patients),
        unique_samples=len(samples),
        errors=tuple(errors),
        warnings=tuple(warnings),
    )
