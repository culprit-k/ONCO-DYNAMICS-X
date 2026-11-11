"""Aggregate-only molecular marker availability audit across breast-cancer cohorts.

This module inventories a small, predeclared descriptive marker panel in public
cBioPortal expression profiles. It does not fit a predictor, harmonize scales,
or export row-level sample/patient identifiers.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from onco_dynamics_x.data_engine.cbioportal import CBioPortalClient


@dataclass(frozen=True)
class CohortProfileSpec:
    """Public study, mRNA molecular profile and compatible sample-list IDs."""

    label: str
    study_id: str
    molecular_profile_id: str
    sample_list_id: str


MARKER_PANEL: tuple[tuple[str, int, str], ...] = (
    ("ESR1", 2099, "Hormone receptor signalling"),
    ("PGR", 5241, "Hormone receptor signalling"),
    ("ERBB2", 2064, "Receptor tyrosine kinase signalling"),
    ("MKI67", 4288, "Proliferation"),
    ("TP53", 7157, "Tumour suppressor"),
    ("BRCA1", 672, "DNA repair"),
    ("BRCA2", 675, "DNA repair"),
    ("EGFR", 1956, "Receptor tyrosine kinase signalling"),
    ("HIF1A", 3091, "Hypoxia response"),
    ("VEGFA", 7422, "Angiogenesis"),
    ("CA9", 768, "Hypoxia response"),
    ("SLC2A1", 6513, "Glucose transport"),
    ("STAT1", 6772, "Interferon signalling"),
    ("CD274", 29126, "Immune checkpoint"),
    ("CXCL9", 4283, "Immune chemokine"),
    ("CD8A", 925, "T-cell marker"),
)

DEFAULT_COHORTS = (
    CohortProfileSpec(
        label="TCGA-BRCA PanCancer Atlas",
        study_id="brca_tcga_pan_can_atlas_2018",
        molecular_profile_id="brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
        sample_list_id="brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
    ),
    CohortProfileSpec(
        label="METABRIC",
        study_id="brca_metabric",
        molecular_profile_id="brca_metabric_mrna",
        sample_list_id="brca_metabric_mrna",
    ),
)


def _summarise_gene_records(
    records: list[dict[str, Any]],
    expected_sample_count: int,
) -> dict[str, Any]:
    """Calculate sample coverage and distribution summaries without IDs."""
    if not isinstance(records, list):
        raise ValueError("Molecular data endpoint must return a JSON list.")

    sample_ids: set[str] = set()
    duplicate_sample_records = 0
    values: list[float] = []
    missing_values = 0
    non_numeric_values = 0
    samples_without_id = 0

    for row in records:
        sample_id = row.get("sampleId")
        if not sample_id:
            samples_without_id += 1
        elif sample_id in sample_ids:
            duplicate_sample_records += 1
        else:
            sample_ids.add(str(sample_id))

        raw_value = row.get("value")
        if raw_value is None:
            missing_values += 1
            continue
        try:
            value = float(raw_value)
        except (TypeError, ValueError):
            non_numeric_values += 1
            continue
        if not np.isfinite(value):
            non_numeric_values += 1
            continue
        values.append(value)

    array = np.asarray(values, dtype=float)
    distribution = None
    if array.size:
        quantiles = np.quantile(array, [0.10, 0.25, 0.50, 0.75, 0.90])
        distribution = {
            "minimum": float(np.min(array)),
            "p10": float(quantiles[0]),
            "p25": float(quantiles[1]),
            "median": float(quantiles[2]),
            "p75": float(quantiles[3]),
            "p90": float(quantiles[4]),
            "maximum": float(np.max(array)),
            "mean": float(np.mean(array)),
            "standard_deviation": float(np.std(array, ddof=0)),
        }

    unique_samples = len(sample_ids)
    return {
        "records_returned": len(records),
        "unique_samples_with_record": unique_samples,
        "sample_coverage_fraction": (
            float(unique_samples / expected_sample_count)
            if expected_sample_count
            else None
        ),
        "missing_value_records": missing_values,
        "non_numeric_or_nonfinite_value_records": non_numeric_values,
        "records_without_sample_identifier": samples_without_id,
        "duplicate_sample_records": duplicate_sample_records,
        "finite_value_count": int(array.size),
        "distribution": distribution,
    }


def audit_molecular_profile(
    client: CBioPortalClient,
    cohort: CohortProfileSpec,
    marker_panel: tuple[tuple[str, int, str], ...] = MARKER_PANEL,
) -> dict[str, Any]:
    """Fetch marker values one gene at a time and emit aggregate statistics only."""
    study = client.get_json(f"studies/{cohort.study_id}")
    profile = client.get_json(
        f"molecular-profiles/{cohort.molecular_profile_id}"
    )
    sample_list = client.get_json(f"sample-lists/{cohort.sample_list_id}")

    if study.get("studyId") != cohort.study_id:
        raise ValueError(f"Unexpected cBioPortal study response for {cohort.label}.")
    if profile.get("molecularProfileId") != cohort.molecular_profile_id:
        raise ValueError(f"Unexpected molecular profile response for {cohort.label}.")
    if sample_list.get("sampleListId") != cohort.sample_list_id:
        raise ValueError(f"Unexpected sample-list response for {cohort.label}.")
    if profile.get("studyId") not in {None, cohort.study_id}:
        raise ValueError(f"Molecular profile does not belong to {cohort.label}.")
    if sample_list.get("studyId") not in {None, cohort.study_id}:
        raise ValueError(f"Sample list does not belong to {cohort.label}.")

    sample_ids = sample_list.get("sampleIds")
    if not isinstance(sample_ids, list) or not sample_ids:
        raise ValueError(f"Sample list is empty or invalid for {cohort.label}.")

    profile_type = str(profile.get("datatype") or "").lower()
    alteration_type = str(profile.get("molecularAlterationType") or "").lower()
    if "mrna" not in profile_type and "mrna" not in alteration_type:
        raise ValueError(
            f"Expected an mRNA expression profile for {cohort.label}; "
            "received a different assay type."
        )

    marker_results = []
    for symbol, entrez_id, biological_context in marker_panel:
        gene = client.get_json(f"genes/{entrez_id}")
        returned_symbol = gene.get("hugoGeneSymbol")
        if returned_symbol and returned_symbol != symbol:
            raise ValueError(
                f"Entrez ID {entrez_id} resolves to {returned_symbol}, not {symbol}."
            )
        if not returned_symbol:
            raise ValueError(f"Entrez ID {entrez_id} has no returned gene symbol.")

        records = client.get_json(
            f"molecular-profiles/{cohort.molecular_profile_id}/molecular-data",
            {
                "sampleListId": cohort.sample_list_id,
                "entrezGeneId": entrez_id,
                "projection": "SUMMARY",
            },
        )
        marker_results.append(
            {
                "gene_symbol": symbol,
                "entrez_gene_id": int(entrez_id),
                "biological_context": biological_context,
                **_summarise_gene_records(records, len(sample_ids)),
            }
        )

    available_markers = sum(
        item["finite_value_count"] > 0 for item in marker_results
    )
    return {
        "cohort_label": cohort.label,
        "study_id": cohort.study_id,
        "molecular_profile_id": cohort.molecular_profile_id,
        "profile_name": profile.get("name"),
        "profile_description": profile.get("description"),
        "data_type": profile.get("datatype"),
        "molecular_alteration_type": profile.get("molecularAlterationType"),
        "reference_genome": study.get("referenceGenome"),
        "sample_list_id": cohort.sample_list_id,
        "sample_list_sample_count": len(sample_ids),
        "marker_panel_size": len(marker_panel),
        "markers_with_finite_values": available_markers,
        "markers_without_finite_values": len(marker_panel) - available_markers,
        "markers": marker_results,
    }


def audit_breast_cancer_marker_panel(
    client: CBioPortalClient | None = None,
    cohorts: tuple[CohortProfileSpec, ...] = DEFAULT_COHORTS,
    marker_panel: tuple[tuple[str, int, str], ...] = MARKER_PANEL,
) -> dict[str, Any]:
    """Audit a fixed marker panel in TCGA and METABRIC without exporting rows."""
    api = client or CBioPortalClient(timeout=60, retries=3)
    if not cohorts:
        raise ValueError("At least one cohort must be supplied.")
    if len({cohort.label for cohort in cohorts}) != len(cohorts):
        raise ValueError("Cohort labels must be unique.")
    if len({symbol for symbol, _, _ in marker_panel}) != len(marker_panel):
        raise ValueError("Marker symbols must be unique.")
    if len({entrez_id for _, entrez_id, _ in marker_panel}) != len(marker_panel):
        raise ValueError("Entrez gene IDs must be unique.")

    cohort_results = [
        audit_molecular_profile(api, cohort, marker_panel)
        for cohort in cohorts
    ]
    by_symbol = {
        result["cohort_label"]: {
            marker["gene_symbol"]: marker
            for marker in result["markers"]
        }
        for result in cohort_results
    }
    common_markers = []
    for symbol, entrez_id, _ in marker_panel:
        if all(
            by_symbol[cohort.label][symbol]["finite_value_count"] > 0
            for cohort in cohorts
        ):
            common_markers.append(
                {"gene_symbol": symbol, "entrez_gene_id": int(entrez_id)}
            )

    return {
        "analysis": "public_breast_cancer_expression_marker_availability_audit",
        "panel": [
            {
                "gene_symbol": symbol,
                "entrez_gene_id": int(entrez_id),
                "biological_context": context,
            }
            for symbol, entrez_id, context in marker_panel
        ],
        "cohorts": cohort_results,
        "cross_cohort_overlap": {
            "marker_count_present_in_all_profiles": len(common_markers),
            "markers_present_in_all_profiles": common_markers,
        },
        "interpretation": {
            "this_is_a_predictive_signature": False,
            "cross_platform_harmonization_performed": False,
            "outcomes_used": False,
            "model_fitted": False,
            "biological_validation_performed": False,
            "patient_or_sample_level_rows_exported": False,
            "raw_expression_values_exported": False,
            "note": (
                "Descriptive marker availability and distribution audit only. "
                "TCGA RNA-seq and METABRIC microarray scales are not directly "
                "comparable. The marker panel is not a validated prognostic, "
                "diagnostic or treatment-response signature."
            ),
        },
        "data_governance": {
            "patient_level_records_written": False,
            "patient_or_sample_identifiers_written": False,
            "raw_expression_values_written": False,
            "aggregate_metrics_only": True,
        },
    }
