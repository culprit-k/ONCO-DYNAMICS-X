"""Public-only NCI Genomic Data Commons metadata query and download client."""
from __future__ import annotations

import csv
import hashlib
import http.client
import json
import os
import re
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any

GDC_API = "https://api.gdc.cancer.gov"
_UUID_PATTERN = re.compile(
    r"^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$"
)


@dataclass(frozen=True)
class GDCFile:
    """Minimal GDC file metadata for provenance and safe open-file download."""

    file_id: str
    file_name: str
    md5sum: str
    file_size: int
    access: str
    data_type: str
    workflow_type: str
    case_ids: tuple[str, ...]
    patient_submitter_ids: tuple[str, ...]
    sample_submitter_ids: tuple[str, ...]
    associated_samples: tuple[tuple[str, str, str], ...]
    sample_types: tuple[str, ...] = ()
    associated_sample_metadata: tuple[tuple[str, str, str, str], ...] = ()


class GDCClient:
    """Small standard-library GDC API client, restricted to open-access queries."""

    def __init__(self, base_url: str = GDC_API, timeout: float = 60.0) -> None:
        if not base_url.startswith("https://"):
            raise ValueError("GDC base URL must use HTTPS.")
        if timeout <= 0:
            raise ValueError("timeout must be positive.")
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout

    def _post_json(self, endpoint: str, payload: dict[str, Any]) -> dict[str, Any]:
        body = json.dumps(payload).encode("utf-8")
        request = urllib.request.Request(
            f"{self.base_url}/{endpoint.lstrip('/')}",
            data=body,
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            method="POST",
        )
        try:
            with urllib.request.urlopen(request, timeout=self.timeout) as response:
                decoded = json.loads(response.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as error:
            raise RuntimeError(f"GDC API request failed: {error}") from error
        if not isinstance(decoded, dict):
            raise RuntimeError("GDC API returned an unexpected response shape.")
        return decoded

    @staticmethod
    def _parse_file(hit: dict[str, Any]) -> GDCFile:
        cases = hit.get("cases") or []
        case_ids: list[str] = []
        patient_ids: list[str] = []
        sample_ids: list[str] = []
        sample_types: list[str] = []
        associated_samples: list[tuple[str, str, str]] = []
        associated_sample_metadata: list[tuple[str, str, str, str]] = []
        for case in cases:
            case_id = str(case.get("case_id") or "")
            patient_id = str(case.get("submitter_id") or "")
            if case_id:
                case_ids.append(case_id)
            if patient_id:
                patient_ids.append(patient_id)
            for sample in case.get("samples") or []:
                sample_id = str(sample.get("submitter_id") or "")
                sample_type = str(sample.get("sample_type") or "")
                if sample_id:
                    sample_ids.append(sample_id)
                    associated_samples.append((case_id, patient_id, sample_id))
                    associated_sample_metadata.append(
                        (case_id, patient_id, sample_id, sample_type)
                    )
                    if sample_type:
                        sample_types.append(sample_type)
        return GDCFile(
            file_id=str(hit.get("file_id") or hit.get("id") or ""),
            file_name=str(hit.get("file_name") or ""),
            md5sum=str(hit.get("md5sum") or ""),
            file_size=int(hit.get("file_size") or 0),
            access=str(hit.get("access") or "").lower(),
            data_type=str(hit.get("data_type") or ""),
            workflow_type=str((hit.get("analysis") or {}).get("workflow_type") or ""),
            case_ids=tuple(case_ids),
            patient_submitter_ids=tuple(patient_ids),
            sample_submitter_ids=tuple(sample_ids),
            associated_samples=tuple(associated_samples),
            sample_types=tuple(sample_types),
            associated_sample_metadata=tuple(associated_sample_metadata),
        )

    def query_open_files(
        self,
        project_id: str = "TCGA-BRCA",
        data_type: str = "Gene Expression Quantification",
        workflow_type: str | None = None,
        page_size: int = 100,
        max_results: int | None = 100,
        sample_type: str | None = None,
    ) -> list[GDCFile]:
        """Search open file metadata, optionally restricting the linked specimen type.

        Results are metadata only. This method does not download files or request
        controlled-access records. Set max_results=None to page through all matches.
        """
        if not project_id.strip() or not data_type.strip():
            raise ValueError("project_id and data_type must be non-empty.")
        if sample_type is not None and not sample_type.strip():
            raise ValueError("sample_type must be non-empty when supplied.")
        if not 1 <= page_size <= 500:
            raise ValueError("page_size must be between 1 and 500.")
        if max_results is not None and max_results < 1:
            raise ValueError("max_results must be positive or None.")

        clauses: list[dict[str, Any]] = [
            {
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": [project_id],
                },
            },
            {
                "op": "=",
                "content": {"field": "files.data_type", "value": data_type},
            },
            {
                "op": "=",
                "content": {"field": "files.access", "value": "open"},
            },
        ]
        if workflow_type:
            clauses.append(
                {
                    "op": "=",
                    "content": {
                        "field": "analysis.workflow_type",
                        "value": workflow_type,
                    },
                }
            )
        if sample_type:
            clauses.append(
                {
                    "op": "=",
                    "content": {
                        "field": "cases.samples.sample_type",
                        "value": sample_type,
                    },
                }
            )

        fields = [
            "file_id",
            "file_name",
            "md5sum",
            "file_size",
            "access",
            "data_type",
            "analysis.workflow_type",
            "cases.case_id",
            "cases.submitter_id",
            "cases.samples.submitter_id",
            "cases.samples.sample_type",
        ]
        results: list[GDCFile] = []
        offset = 0
        total = None
        while total is None or offset < total:
            remaining = (
                page_size
                if max_results is None
                else min(page_size, max_results - len(results))
            )
            if remaining <= 0:
                break
            payload = {
                "filters": {"op": "and", "content": clauses},
                "format": "JSON",
                "fields": ",".join(fields),
                "size": remaining,
                "from": offset,
            }
            response = self._post_json("files", payload)
            data = response.get("data") or {}
            hits = data.get("hits") or []
            pagination = data.get("pagination") or {}
            if total is None:
                total = int(pagination.get("total", len(hits)))
            parsed = [self._parse_file(hit) for hit in hits]
            for item in parsed:
                if item.access != "open":
                    raise RuntimeError(
                        "GDC response included a non-open file despite the open filter."
                    )
                if not item.file_id or not item.file_name:
                    raise RuntimeError("GDC response omitted file identity metadata.")
                if sample_type and (
                    not item.sample_types
                    or any(value != sample_type for value in item.sample_types)
                ):
                    raise RuntimeError(
                        "GDC response included a file without the requested linked sample type."
                    )
            results.extend(parsed)
            if not hits:
                break
            offset += len(hits)
        return results

    def download_open_file(
        self,
        metadata: GDCFile,
        destination: str | Path,
        retries: int = 4,
        retry_delay: float = 1.0,
    ) -> Path:
        """Download one open file, retry transient failures, and verify size and MD5."""
        if metadata.access != "open":
            raise PermissionError("Only open-access GDC files can be downloaded by this client.")
        if not _UUID_PATTERN.fullmatch(metadata.file_id):
            raise ValueError("file_id must be a valid GDC UUID.")
        if not metadata.file_name or Path(metadata.file_name).name != metadata.file_name:
            raise ValueError("file_name must be a plain filename without path components.")
        if metadata.file_size < 0 or not re.fullmatch(r"[0-9a-fA-F]{32}", metadata.md5sum):
            raise ValueError("Expected file size and MD5 checksum are required.")
        if retries < 0:
            raise ValueError("retries must be non-negative.")
        if retry_delay < 0:
            raise ValueError("retry_delay must be non-negative.")

        destination_path = Path(destination)
        destination_path.mkdir(parents=True, exist_ok=True)
        final_path = destination_path / metadata.file_name
        partial_path = final_path.with_suffix(final_path.suffix + ".partial")
        request = urllib.request.Request(
            f"{self.base_url}/data/{metadata.file_id}",
            headers={
                "Accept": "application/octet-stream",
                "User-Agent": "ONCO-DYNAMICS-X/0.1 research data client",
            },
        )
        last_error: Exception | None = None
        for attempt in range(retries + 1):
            digest = hashlib.md5(usedforsecurity=False)
            bytes_written = 0
            try:
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    with partial_path.open("wb") as output:
                        while chunk := response.read(1024 * 1024):
                            output.write(chunk)
                            digest.update(chunk)
                            bytes_written += len(chunk)
                if bytes_written != metadata.file_size:
                    raise RuntimeError(
                        f"Downloaded size {bytes_written} does not match metadata "
                        f"size {metadata.file_size}."
                    )
                if digest.hexdigest().lower() != metadata.md5sum.lower():
                    raise RuntimeError("Downloaded file MD5 does not match GDC metadata.")
                os.replace(partial_path, final_path)
                return final_path
            except (
                OSError,
                urllib.error.URLError,
                http.client.HTTPException,
                TimeoutError,
                RuntimeError,
            ) as error:
                last_error = error
                partial_path.unlink(missing_ok=True)
                if attempt >= retries:
                    break
                if retry_delay:
                    time.sleep(retry_delay * (2**attempt))
        raise RuntimeError(
            f"Could not download and verify GDC file {metadata.file_id} "
            f"after {retries + 1} attempts: {last_error}"
        ) from last_error


def write_sample_manifest(
    files: list[GDCFile],
    output: str | Path,
    cohort_id: str,
    modality: str,
) -> Path:
    """Write a sample-level manifest with explicit case-to-sample linkages."""
    if not cohort_id.strip() or not modality.strip():
        raise ValueError("cohort_id and modality must be non-empty.")
    manifest_path = Path(output)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "cohort_id",
        "patient_id",
        "specimen_id",
        "sample_id",
        "sample_type",
        "assay_type",
        "modality",
        "data_access",
        "source_accession",
    ]
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        seen: set[tuple[str, str, str]] = set()
        for file in files:
            if file.access != "open":
                raise PermissionError("Manifest export accepts open-access file metadata only.")
            sample_type_by_id = {
                sample_id: sample_type
                for _, _, sample_id, sample_type in file.associated_sample_metadata
            }
            for case_id, patient_id, sample_id in file.associated_samples:
                association = (case_id, patient_id, sample_id)
                if association in seen:
                    continue
                seen.add(association)
                writer.writerow(
                    {
                        "cohort_id": cohort_id,
                        "patient_id": patient_id,
                        "specimen_id": sample_id,
                        "sample_id": sample_id,
                        "sample_type": sample_type_by_id.get(sample_id, ""),
                        "assay_type": file.data_type,
                        "modality": modality,
                        "data_access": file.access,
                        "source_accession": file.file_id,
                    }
                )
    return manifest_path
