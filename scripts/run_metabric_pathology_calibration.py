"""Evaluate fixed molecular modules against METABRIC pathology annotations.

Only aggregate metrics are written. The cohort is cross-sectional; these
associations do not calibrate spatial PDE rates or establish causal biology.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from onco_dynamics_x.calibration.pathology import (
    evaluate_module_pathology_association,
)
from onco_dynamics_x.data_engine.cbioportal import CBioPortalClient
from onco_dynamics_x.data_engine.cohort_metadata import (
    parse_cbioportal_sample_metadata,
)
from onco_dynamics_x.data_engine.molecular_mechanistic_bridge import (
    score_molecular_modules,
)
from onco_dynamics_x.data_engine import molecular_model
from onco_dynamics_x.data_engine.molecular_model import (
    PAM50_GENES,
    rank_normalise_expression,
)

STUDY_ID = "brca_metabric"
PROFILE_ID = "brca_metabric_mrna"
SAMPLE_LIST_ID = "brca_metabric_mrna"
TARGETS = {
    "ER_STATUS": "SAMPLE",
    "PR_STATUS": "SAMPLE",
    "HER2_STATUS": "SAMPLE",
    "GRADE": "SAMPLE",
    "CELLULARITY": "PATIENT",
}
CHUNK_SIZE = 300


def _resolve_entrez_ids(
    client: CBioPortalClient,
) -> dict[str, int]:
    """Resolve the fixed panel through current cBioPortal gene records."""
    resolved: dict[str, int] = {}
    for symbol in PAM50_GENES:
        last_error: Exception | None = None
        for candidate in molecular_model.PAM50_SYMBOL_ALIASES.get(symbol, (symbol,)):
            try:
                gene = client.get_json(f"genes/{candidate}")
            except Exception as error:
                last_error = error
                continue
            entrez_id = gene.get("entrezGeneId")
            if entrez_id:
                resolved[symbol] = int(entrez_id)
                break
        if symbol not in resolved:
            raise RuntimeError(f"Could not resolve {symbol}: {last_error}")
    if len(set(resolved.values())) != len(resolved):
        raise ValueError("Panel symbols unexpectedly resolve to duplicate Entrez IDs.")
    return resolved


def _fetch_expression_records(
    client: CBioPortalClient,
    sample_ids: list[str],
    entrez_ids: list[int],
) -> tuple[np.ndarray, dict[str, int]]:
    """Fetch fixed-panel expression in bounded batches."""
    sample_index = {sample_id: index for index, sample_id in enumerate(sample_ids)}
    gene_index = {gene_id: index for index, gene_id in enumerate(entrez_ids)}
    values = np.full((len(sample_ids), len(entrez_ids)), np.nan, dtype=float)
    record_count = 0
    duplicates = 0
    for start in range(0, len(sample_ids), CHUNK_SIZE):
        records = client.post_json(
            f"molecular-profiles/{PROFILE_ID}/molecular-data/fetch",
            {
                "sampleIds": sample_ids[start : start + CHUNK_SIZE],
                "entrezGeneIds": entrez_ids,
            },
        )
        if not isinstance(records, list):
            raise RuntimeError("cBioPortal molecular-data/fetch returned a non-list.")
        seen: set[tuple[int, int]] = set()
        for record in records:
            record_count += 1
            sample_id = str(record.get("sampleId") or "")
            try:
                gene_id = int(record.get("entrezGeneId"))
                value = float(record.get("value"))
            except (TypeError, ValueError):
                continue
            if sample_id not in sample_index or gene_id not in gene_index:
                continue
            pair = (sample_index[sample_id], gene_index[gene_id])
            if pair in seen:
                duplicates += 1
                values[pair] = np.nan
                continue
            seen.add(pair)
            if np.isfinite(value):
                values[pair] = value
    qc = {
        "records_returned": int(record_count),
        "duplicate_records": int(duplicates),
        "finite_measurements": int(np.isfinite(values).sum()),
        "expected_measurements": int(values.size),
    }
    return values, qc




def _unique_patient_annotations(
    records: list[dict[str, Any]],
    attribute_id: str,
    level: str,
    sample_to_patient: dict[str, str],
) -> tuple[dict[str, str], int]:
    grouped: dict[str, set[str]] = defaultdict(set)
    entity_key = "patientId" if level == "PATIENT" else "sampleId"
    for row in records:
        value = row.get("value")
        if row.get("clinicalAttributeId") != attribute_id:
            continue
        if value is None or not str(value).strip():
            continue
        entity_id = str(row.get(entity_key) or "")
        patient_id = (
            entity_id
            if level == "PATIENT"
            else sample_to_patient.get(entity_id, "")
        )
        if patient_id:
            grouped[patient_id].add(str(value).strip())
    values = {
        patient_id: next(iter(observed))
        for patient_id, observed in grouped.items()
        if len(observed) == 1
    }
    ambiguous = sum(len(observed) > 1 for observed in grouped.values())
    return values, int(ambiguous)


def run_metabric_pathology_calibration(
    *,
    n_splits: int = 5,
    n_bootstrap: int = 500,
    seed: int = 42,
) -> dict[str, Any]:
    """Fetch METABRIC primary-tumour expression and evaluate module associations."""
    client = CBioPortalClient(timeout=60, retries=3)
    gene_ids = _resolve_entrez_ids(client)
    sample_list = client.get_json(f"sample-lists/{SAMPLE_LIST_ID}")
    sample_ids = sample_list.get("sampleIds")
    if not isinstance(sample_ids, list) or not sample_ids:
        raise RuntimeError("METABRIC expression sample list is empty or invalid.")

    samples = client.get_all(
        f"studies/{STUDY_ID}/samples", {"projection": "DETAILED"}
    )
    metadata = parse_cbioportal_sample_metadata(samples)
    primary_sample_ids = [
        str(sample_id)
        for sample_id in sample_ids
        if metadata.get(str(sample_id), {}).get("type") == "Primary Solid Tumor"
        and metadata.get(str(sample_id), {}).get("patientId")
    ]
    if not primary_sample_ids:
        raise RuntimeError("No verified Primary Solid Tumor samples were found.")
    sample_to_patient = {
        sample_id: str(metadata[sample_id]["patientId"])
        for sample_id in primary_sample_ids
    }

    raw_expression, expression_qc = _fetch_expression_records(
        client, primary_sample_ids, [gene_ids[gene] for gene in PAM50_GENES]
    )
    patient_rows: dict[str, list[np.ndarray]] = defaultdict(list)
    for index, sample_id in enumerate(primary_sample_ids):
        patient_rows[sample_to_patient[sample_id]].append(raw_expression[index])
    patient_ids = sorted(patient_rows)
    patient_expression = np.stack(
        [np.mean(patient_rows[patient_id], axis=0) for patient_id in patient_ids]
    )
    rank_result = rank_normalise_expression(
        patient_expression,
        PAM50_GENES,
        minimum_observed_fraction=0.95,
    )
    module_profile = score_molecular_modules(
        rank_result.values, PAM50_GENES, minimum_gene_coverage=0.75
    )

    patient_records = client.get_all(
        f"studies/{STUDY_ID}/clinical-data",
        {"clinicalDataType": "PATIENT"},
    )
    sample_records = client.get_all(
        f"studies/{STUDY_ID}/clinical-data",
        {"clinicalDataType": "SAMPLE"},
    )
    reports: dict[str, Any] = {}
    for attribute_id, level in TARGETS.items():
        records = patient_records if level == "PATIENT" else sample_records
        labels, ambiguous = _unique_patient_annotations(
            records, attribute_id, level, sample_to_patient
        )
        report = evaluate_module_pathology_association(
            module_profile.sample_scores,
            patient_ids,
            labels,
            module_names=module_profile.module_names,
            n_splits=n_splits,
            n_bootstrap=n_bootstrap,
            seed=seed,
        )
        report["source_attribute"] = attribute_id
        report["source_level"] = level
        report["ambiguous_patients_excluded"] = ambiguous
        reports[attribute_id] = report

    return {
        "analysis": "METABRIC_cross_validated_molecular_module_pathology_association",
        "cohort": {
            "study_id": STUDY_ID,
            "molecular_profile_id": PROFILE_ID,
            "expression_sample_list_count": len(sample_ids),
            "verified_primary_solid_tumour_samples": len(primary_sample_ids),
            "unique_primary_tumour_patients": len(patient_ids),
            "patients_with_multiple_primary_samples": sum(
                len(rows) > 1 for rows in patient_rows.values()
            ),
            "expression_fetch_qc": expression_qc,
        },
        "feature_construction": {
            "panel": "fixed 50-gene PAM50 research inventory",
            "transformation": "within-cohort mid-rank inverse-normal per gene",
            "modules": list(module_profile.module_names),
            "module_gene_coverage": module_profile.gene_coverage,
            "uses_pathology_labels_for_features": False,
            "transductive_expression_transform": True,
        },
        "pathology_targets": reports,
        "limitations": [
            "Cross-sectional association with existing METABRIC clinical annotations only.",
            "No Ki-67 field was identified in the aggregate target audit.",
            "Field definitions and assay provenance require source-publication review.",
            "Cohort-wise rank normalisation is transductive, not prospective single-sample processing.",
            "This analysis does not estimate PDE parameters or establish causal mechanisms.",
            "No patient identifiers, row-level annotations, or expression values are exported.",
        ],
        "data_governance": {
            "patient_ids_exported": False,
            "sample_ids_exported": False,
            "row_level_annotations_exported": False,
            "raw_expression_exported": False,
            "aggregate_metrics_only": True,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="results/biological_calibration/metabric_pathology_module_association.json",
    )
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if args.folds < 2 or args.bootstrap < 0:
        parser.error("--folds must be >=2 and --bootstrap must be non-negative")
    report = run_metabric_pathology_calibration(
        n_splits=args.folds,
        n_bootstrap=args.bootstrap,
        seed=args.seed,
    )
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    for target, result in report["pathology_targets"].items():
        metrics = result["metrics"]
        print(
            target,
            "n=", result["patients_evaluated"],
            "balanced_accuracy=", round(metrics["balanced_accuracy"], 4),
            "macro_f1=", round(metrics["macro_f1"], 4),
        )
    print(f"Aggregate report: {path}")


if __name__ == "__main__":
    main()
