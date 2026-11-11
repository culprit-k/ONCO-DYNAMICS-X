"""Train a fixed-panel TCGA molecular model and evaluate it once in METABRIC.

The model uses the published PAM50 50-gene inventory as a fixed research
feature panel, but does not copy PAM50 centroids or claim to reproduce the
PAM50/Prosigna classifier. Expression is rank-normalised within each cohort
without using clinical outcomes, then TCGA-fitted PCA loadings are frozen and
applied to METABRIC. This is an exploratory, cohort-level external transport
experiment; it is not suitable for individual-patient deployment.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from onco_dynamics_x.data_engine import molecular_model
from onco_dynamics_x.data_engine.cbioportal import CBioPortalClient
from onco_dynamics_x.data_engine.cohort_metadata import parse_cbioportal_sample_metadata
from onco_dynamics_x.data_engine.survival_analysis import fetch_survival_cohort
from onco_dynamics_x.data_engine.survival_metrics import evaluate_ipcw_survival

COHORTS = (
    {
        "label": "TCGA-BRCA PanCancer Atlas",
        "study_id": "brca_tcga_pan_can_atlas_2018",
        "profile_id": "brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
        "sample_list_id": "brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
        "age_attribute": "AGE",
        "stage_attribute": "AJCC_PATHOLOGIC_TUMOR_STAGE",
        "stage_level": "PATIENT",
    },
    {
        "label": "METABRIC",
        "study_id": "brca_metabric",
        "profile_id": "brca_metabric_mrna",
        "sample_list_id": "brca_metabric_mrna",
        "age_attribute": "AGE_AT_DIAGNOSIS",
        "stage_attribute": "TUMOR_STAGE",
        "stage_level": "SAMPLE",
    },
)

CHUNK_SIZE = 300
N_COMPONENTS = 5
COX_PENALIZER = 0.5
HORIZONS_MONTHS = (12.0, 36.0, 60.0)


def _resolve_entrez_ids(
    client: CBioPortalClient,
) -> tuple[dict[str, int], dict[str, str]]:
    """Resolve panel symbols through cBioPortal, allowing historical aliases."""
    resolved_ids: dict[str, int] = {}
    resolved_symbols: dict[str, str] = {}
    for panel_symbol in molecular_model.PAM50_GENES:
        candidates = molecular_model.PAM50_SYMBOL_ALIASES.get(panel_symbol, (panel_symbol,))
        last_error: Exception | None = None
        for candidate in candidates:
            try:
                gene = client.get_json(f"genes/{candidate}")
            except Exception as error:  # API returns 404 for retired HUGO aliases.
                last_error = error
                continue
            entrez_id = gene.get("entrezGeneId")
            returned_symbol = gene.get("hugoGeneSymbol")
            if not entrez_id or not returned_symbol:
                continue
            resolved_ids[panel_symbol] = int(entrez_id)
            resolved_symbols[panel_symbol] = str(returned_symbol)
            break
        if panel_symbol not in resolved_ids:
            raise RuntimeError(
                f"Could not resolve PAM50 panel gene {panel_symbol}: {last_error}"
            )
    if len(set(resolved_ids.values())) != len(resolved_ids):
        raise ValueError("Two PAM50 panel symbols unexpectedly resolved to one Entrez ID.")
    return resolved_ids, resolved_symbols


def _is_primary_tumour(sample_type: Any) -> bool:
    token = " ".join(str(sample_type or "").lower().replace("_", " ").split())
    return "primary" in token and "tumor" in token


def _fetch_expression_records(
    client: CBioPortalClient,
    profile_id: str,
    sample_ids: list[str],
    entrez_ids: list[int],
) -> tuple[np.ndarray, dict[str, Any]]:
    """Fetch the gene-by-sample matrix in bounded sample batches."""
    sample_index = {sample_id: index for index, sample_id in enumerate(sample_ids)}
    entrez_index = {entrez_id: index for index, entrez_id in enumerate(entrez_ids)}
    values = np.full((len(sample_ids), len(entrez_ids)), np.nan, dtype=float)
    duplicate_records = 0
    unrecognised_records = 0
    nonfinite_records = 0
    record_count = 0
    seen_sample_gene_pairs: set[tuple[int, int]] = set()

    for start in range(0, len(sample_ids), CHUNK_SIZE):
        sample_chunk = sample_ids[start : start + CHUNK_SIZE]
        records = client.post_json(
            f"molecular-profiles/{profile_id}/molecular-data/fetch",
            {
                "sampleIds": sample_chunk,
                "entrezGeneIds": entrez_ids,
            },
        )
        if not isinstance(records, list):
            raise RuntimeError("cBioPortal molecular-data/fetch returned a non-list response.")
        for record in records:
            record_count += 1
            sample_id = str(record.get("sampleId") or "")
            try:
                entrez_id = int(record.get("entrezGeneId"))
            except (TypeError, ValueError):
                unrecognised_records += 1
                continue
            if sample_id not in sample_index or entrez_id not in entrez_index:
                unrecognised_records += 1
                continue
            row = sample_index[sample_id]
            column = entrez_index[entrez_id]
            pair = (row, column)
            if pair in seen_sample_gene_pairs:
                duplicate_records += 1
                values[row, column] = np.nan
                continue
            seen_sample_gene_pairs.add(pair)
            raw_value = record.get("value")
            try:
                numeric_value = float(raw_value)
            except (TypeError, ValueError):
                nonfinite_records += 1
                continue
            if not np.isfinite(numeric_value):
                nonfinite_records += 1
                continue
            values[row, column] = numeric_value

    return values, {
        "records_returned": int(record_count),
        "duplicate_sample_gene_records": int(duplicate_records),
        "unrecognised_records": int(unrecognised_records),
        "nonfinite_or_unparseable_records": int(nonfinite_records),
        "finite_sample_gene_measurements": int(np.isfinite(values).sum()),
        "expected_sample_gene_measurements": int(values.size),
    }


def _audit_subtype_annotations(
    client: CBioPortalClient,
    study_id: str,
) -> dict[str, Any]:
    """Summarise available molecular-subtype clinical fields without exporting rows."""
    attributes = client.get_all(
        f"studies/{study_id}/clinical-attributes",
        {"projection": "SUMMARY"},
    )
    patient_records = client.get_all(
        f"studies/{study_id}/clinical-data",
        {"clinicalDataType": "PATIENT", "projection": "SUMMARY"},
    )
    sample_records = client.get_all(
        f"studies/{study_id}/clinical-data",
        {"clinicalDataType": "SAMPLE", "projection": "SUMMARY"},
    )
    result = []
    selected_patient_labels: dict[str, str] = {}
    preferred_attribute = (
        "SUBTYPE"
        if study_id == "brca_tcga_pan_can_atlas_2018"
        else "CLAUDIN_SUBTYPE"
        if study_id == "brca_metabric"
        else ""
    )
    for attribute in attributes:
        attribute_id = str(attribute.get("clinicalAttributeId") or "")
        display = str(attribute.get("displayName") or "")
        searchable = f"{attribute_id} {display}".upper()
        if not any(term in searchable for term in ("PAM50", "SUBTYPE", "MOLECULAR_CLASS")):
            continue
        patient_level = bool(attribute.get("patientAttribute"))
        records = patient_records if patient_level else sample_records
        field_records = [
            item for item in records
            if item.get("clinicalAttributeId") == attribute_id
            and item.get("value") is not None
            and str(item.get("value")).strip()
        ]
        counts = Counter(str(item.get("value")).strip() for item in field_records)
        entity_key = "patientId" if patient_level else "sampleId"
        result.append(
            {
                "clinical_attribute_id": attribute_id,
                "display_name": display or attribute_id,
                "level": "PATIENT" if patient_level else "SAMPLE",
                "records_with_value": int(len(field_records)),
                "unique_entities_with_value": int(
                    len({item.get(entity_key) for item in field_records if item.get(entity_key)})
                ),
                "value_counts": dict(sorted(counts.items())),
            }
        )
        if patient_level and attribute_id.upper() == preferred_attribute:
            selected_patient_labels = {
                str(item.get("patientId")): str(item.get("value")).strip()
                for item in field_records
                if item.get("patientId") and str(item.get("value")).strip()
            }
    return {
        "candidate_subtype_fields_found": int(len(result)),
        "fields": result,
        "patient_subtype_labels": selected_patient_labels,
        "selected_subtype_attribute": preferred_attribute or None,
        "note": (
            "These are existing cohort annotations. The selected subtype labels "
            "are used only as supervised training targets in TCGA and held-out "
            "evaluation labels in METABRIC; they are never used to construct "
            "expression features or tune the model."
        ),
    }


def _load_cohort(
    client: CBioPortalClient,
    cohort: dict[str, str],
    gene_ids: dict[str, int],
) -> dict[str, Any]:
    """Load one cohort in memory, align primary-tumour expression to OS records."""
    survival_rows, survival_qc = fetch_survival_cohort(
        client=client,
        study_id=cohort["study_id"],
        sample_list_id=cohort["sample_list_id"],
        age_attribute=cohort["age_attribute"],
        stage_attribute=cohort["stage_attribute"],
        stage_level=cohort["stage_level"],
        require_stage=False,
    )
    survival_by_patient = {row["patient_id"]: row for row in survival_rows}
    samples = client.get_all(
        f"studies/{cohort['study_id']}/samples",
        {"projection": "DETAILED"},
    )
    sample_list = client.get_json(f"sample-lists/{cohort['sample_list_id']}")
    profile = client.get_json(f"molecular-profiles/{cohort['profile_id']}")
    if profile.get("molecularProfileId") != cohort["profile_id"]:
        raise RuntimeError(f"Unexpected molecular profile for {cohort['label']}.")
    if profile.get("studyId") not in {None, cohort["study_id"]}:
        raise RuntimeError(f"Profile does not belong to {cohort['label']}.")
    profile_type = str(profile.get("datatype") or "").lower()
    alteration_type = str(profile.get("molecularAlterationType") or "").lower()
    if "mrna" not in profile_type and "mrna" not in alteration_type:
        raise RuntimeError(
            f"Expected an mRNA expression profile for {cohort['label']}."
        )
    if sample_list.get("sampleListId") != cohort["sample_list_id"]:
        raise RuntimeError(f"Unexpected sample list for {cohort['label']}.")
    sample_ids = sample_list.get("sampleIds")
    if not isinstance(sample_ids, list) or not sample_ids:
        raise RuntimeError(f"Empty or invalid expression sample list for {cohort['label']}.")

    sample_metadata = parse_cbioportal_sample_metadata(samples)
    sample_type_counts = Counter(
        sample_metadata.get(sample_id, {}).get("type") or "UNKNOWN"
        for sample_id in sample_ids
    )
    primary_ids = [
        sample_id
        for sample_id in sample_ids
        if _is_primary_tumour(sample_metadata.get(sample_id, {}).get("type"))
    ]
    unknown_ids = [
        sample_id
        for sample_id in sample_ids
        if not sample_metadata.get(sample_id, {}).get("type")
    ]
    # If the API does not expose sample type for this study, retain the sample
    # list but make the lack of type filtering explicit in QC.
    if not primary_ids and unknown_ids:
        primary_ids = list(unknown_ids)
        sample_type_filter = "sample type unavailable; all sample-list entries retained"
    else:
        sample_type_filter = "Primary Solid Tumor only"

    if not primary_ids:
        raise RuntimeError(
            f"No primary tumour expression samples were identified for {cohort['label']}."
        )
    expression_matrix, expression_qc = _fetch_expression_records(
        client,
        cohort["profile_id"],
        primary_ids,
        [gene_ids[gene] for gene in molecular_model.PAM50_GENES],
    )

    selected_sample_ids: list[str] = []
    selected_patient_ids: list[str] = []
    selected_patient_rows: list[dict[str, Any]] = []
    selected_sample_ids_seen: set[str] = set()
    patients_seen: set[str] = set()
    duplicate_samples_per_patient = 0
    samples_without_patient = 0
    for sample_id in primary_ids:
        metadata = sample_metadata.get(sample_id)
        patient_id = metadata.get("patientId") if metadata else None
        if not patient_id:
            samples_without_patient += 1
            continue
        if patient_id not in survival_by_patient:
            continue
        if patient_id in patients_seen:
            duplicate_samples_per_patient += 1
            continue
        if sample_id in selected_sample_ids_seen:
            continue
        patients_seen.add(patient_id)
        selected_sample_ids_seen.add(sample_id)
        selected_sample_ids.append(sample_id)
        selected_patient_ids.append(str(patient_id))
        selected_patient_rows.append(survival_by_patient[patient_id])

    row_index = {sample_id: index for index, sample_id in enumerate(primary_ids)}
    selected_values = np.asarray(
        [expression_matrix[row_index[sample_id]] for sample_id in selected_sample_ids],
        dtype=float,
    )
    if selected_values.ndim != 2 or selected_values.shape[1] != len(molecular_model.PAM50_GENES):
        raise RuntimeError("Selected patient expression matrix has an invalid shape.")

    subtype_audit = _audit_subtype_annotations(client, cohort["study_id"])
    patient_subtype_labels = subtype_audit.pop("patient_subtype_labels", {})

    return {
        "label": cohort["label"],
        "study_id": cohort["study_id"],
        "profile_id": cohort["profile_id"],
        "profile_name": profile.get("name"),
        "profile_description": profile.get("description"),
        "n_expression_samples_in_list": int(len(sample_ids)),
        "n_primary_tumour_samples": (
            int(len(primary_ids))
            if sample_type_filter == "Primary Solid Tumor only"
            else None
        ),
        "n_expression_samples_used": int(len(primary_ids)),
        "n_unique_primary_tumour_patients_with_os": int(len(selected_patient_rows)),
        "n_os_events": int(sum(int(row["event"]) for row in selected_patient_rows)),
        "n_censored": int(sum(1 - int(row["event"]) for row in selected_patient_rows)),
        "sample_type_counts": dict(sorted(sample_type_counts.items())),
        "sample_type_filter": sample_type_filter,
        "duplicate_primary_samples_excluded": int(duplicate_samples_per_patient),
        "samples_without_patient_link": int(samples_without_patient),
        "survival_assembly_qc": survival_qc,
        "expression_fetch_qc": expression_qc,
        "subtype_annotation_audit": subtype_audit,
        "patient_subtype_labels": patient_subtype_labels,
        "patient_ids": selected_patient_ids,
        "patient_rows": selected_patient_rows,
        "sample_ids": selected_sample_ids,
        "expression_values": selected_values,
        "raw_expression_missing_fraction_by_gene": {
            gene: float(1.0 - np.isfinite(selected_values[:, index]).mean())
            for index, gene in enumerate(molecular_model.PAM50_GENES)
        },
    }


def _fit_cox_model(
    train_rows: list[dict[str, Any]],
    test_rows: list[dict[str, Any]],
    train_pcs: np.ndarray | None,
    test_pcs: np.ndarray | None,
    *,
    penalizer: float,
    n_bootstrap: int,
    seed: int,
) -> dict[str, Any]:
    from lifelines import CoxPHFitter
    from lifelines.statistics import proportional_hazard_test
    from lifelines.utils import concordance_index

    train = pd.DataFrame(
        {
            "age_years": [row["age_years"] for row in train_rows],
            "duration_months": [row["duration_months"] for row in train_rows],
            "event": [row["event"] for row in train_rows],
        }
    )
    test = pd.DataFrame(
        {
            "age_years": [row["age_years"] for row in test_rows],
            "duration_months": [row["duration_months"] for row in test_rows],
            "event": [row["event"] for row in test_rows],
        }
    )
    age_mean = float(train["age_years"].mean())
    age_sd = float(train["age_years"].std(ddof=0))
    if not np.isfinite(age_sd) or age_sd <= 0:
        raise ValueError("Training-cohort age variation is invalid.")
    train["age_z"] = (train["age_years"] - age_mean) / age_sd
    test["age_z"] = (test["age_years"] - age_mean) / age_sd

    age_model = CoxPHFitter(penalizer=penalizer)
    age_model.fit(
        train[["age_z", "duration_months", "event"]],
        duration_col="duration_months",
        event_col="event",
    )
    molecular_train = train[["age_z"]].copy()
    molecular_test = test[["age_z"]].copy()
    if train_pcs is not None and test_pcs is not None:
        if train_pcs.shape[0] != len(train) or test_pcs.shape[0] != len(test):
            raise ValueError("PCA scores do not align with survival rows.")
        for index in range(train_pcs.shape[1]):
            feature_name = f"PC{index + 1}"
            molecular_train[feature_name] = train_pcs[:, index]
            molecular_test[feature_name] = test_pcs[:, index]

    molecular_train["duration_months"] = train["duration_months"].to_numpy()
    molecular_train["event"] = train["event"].to_numpy()
    molecular_model = CoxPHFitter(penalizer=penalizer)
    molecular_model.fit(
        molecular_train,
        duration_col="duration_months",
        event_col="event",
    )

    test_durations = test["duration_months"].to_numpy(dtype=float)
    test_events = test["event"].to_numpy(dtype=int)
    age_risk = age_model.predict_partial_hazard(test[["age_z"]]).to_numpy().reshape(-1)
    molecular_risk = (
        molecular_model.predict_partial_hazard(molecular_test).to_numpy().reshape(-1)
    )
    age_c_index = float(concordance_index(test_durations, -age_risk, test_events))
    molecular_c_index = float(
        concordance_index(test_durations, -molecular_risk, test_events)
    )

    rng = np.random.default_rng(seed)
    age_bootstrap: list[float] = []
    molecular_bootstrap: list[float] = []
    delta_bootstrap: list[float] = []
    for _ in range(n_bootstrap):
        indices = rng.integers(0, len(test_durations), size=len(test_durations))
        try:
            age_value = float(
                concordance_index(test_durations[indices], -age_risk[indices], test_events[indices])
            )
            molecular_value = float(
                concordance_index(
                    test_durations[indices], -molecular_risk[indices], test_events[indices]
                )
            )
        except (ZeroDivisionError, ValueError):
            continue
        if np.isfinite(age_value) and np.isfinite(molecular_value):
            age_bootstrap.append(age_value)
            molecular_bootstrap.append(molecular_value)
            delta_bootstrap.append(molecular_value - age_value)

    def ci(values: list[float]) -> list[float] | None:
        if len(values) < max(30, int(0.5 * n_bootstrap)):
            return None
        return [
            float(np.quantile(values, 0.025)),
            float(np.quantile(values, 0.975)),
        ]

    horizon_metrics: dict[str, Any] = {}
    age_survival = age_model.predict_survival_function(
        test[["age_z"]], times=list(HORIZONS_MONTHS)
    ).to_numpy(dtype=float).T
    molecular_survival = molecular_model.predict_survival_function(
        molecular_test, times=list(HORIZONS_MONTHS)
    ).to_numpy(dtype=float).T
    for horizon_index, horizon in enumerate(HORIZONS_MONTHS):
        age_report = evaluate_ipcw_survival(
            test_durations,
            test_events,
            age_survival[:, horizon_index],
            horizons_months=(horizon,),
            n_bootstrap=n_bootstrap,
            seed=seed + horizon_index,
        )["horizons"][0]
        molecular_report = evaluate_ipcw_survival(
            test_durations,
            test_events,
            molecular_survival[:, horizon_index],
            horizons_months=(horizon,),
            n_bootstrap=n_bootstrap,
            seed=seed + 100 + horizon_index,
        )["horizons"][0]
        horizon_metrics[f"{int(horizon)}_months"] = {
            "age_only": age_report,
            "age_plus_pam50_pc": molecular_report,
        }

    ph_test = proportional_hazard_test(
        molecular_model,
        molecular_train,
        time_transform="rank",
    ).summary
    coefficients = []
    for feature, coefficient in molecular_model.summary.iterrows():
        coefficients.append(
            {
                "feature": str(feature),
                "coefficient": float(coefficient["coef"]),
                "hazard_ratio": float(coefficient["exp(coef)"]),
                "hazard_ratio_ci_95": [
                    float(coefficient["exp(coef) lower 95%"]),
                    float(coefficient["exp(coef) upper 95%"]),
                ],
                "p_value": float(coefficient["p"]),
            }
        )
    return {
        "model": "Ridge-penalised Cox model using age and TCGA-fitted PAM50-panel PCs",
        "age_only_reference": {
            "harrell_c_index": age_c_index,
            "harrell_c_index_bootstrap_95ci": ci(age_bootstrap),
        },
        "molecular_model": {
            "harrell_c_index": molecular_c_index,
            "harrell_c_index_bootstrap_95ci": ci(molecular_bootstrap),
            "delta_c_index_vs_age_only": float(molecular_c_index - age_c_index),
            "paired_bootstrap_delta_c_index_95ci": ci(delta_bootstrap),
            "bootstrap_replicates_requested": int(n_bootstrap),
            "bootstrap_replicates_valid": len(delta_bootstrap),
        },
        "time_dependent_metrics": horizon_metrics,
        "model_coefficients": coefficients,
        "proportional_hazards_rank_test_p_values": {
            str(index): float(row["p"])
            for index, row in ph_test.iterrows()
        },
        "training_transform": {
            "age_mean_tcga_years": age_mean,
            "age_sd_tcga_years": age_sd,
            "cox_l2_penalizer": float(penalizer),
            "pca_components": 0 if train_pcs is None else int(train_pcs.shape[1]),
        },
    }


def run_pam50_molecular_validation(
    *,
    n_components: int = N_COMPONENTS,
    n_bootstrap: int = 250,
    seed: int = 42,
    cox_penalizer: float = COX_PENALIZER,
) -> dict[str, Any]:
    """Fetch, harmonise and externally evaluate the fixed 50-gene molecular panel."""
    client = CBioPortalClient(timeout=60, retries=3)
    gene_ids, resolved_symbols = _resolve_entrez_ids(client)
    cohorts = [_load_cohort(client, cohort, gene_ids) for cohort in COHORTS]

    tcga, metabric = cohorts
    common_genes = [
        gene for gene in molecular_model.PAM50_GENES
        if tcga["raw_expression_missing_fraction_by_gene"][gene] <= 0.05
        and metabric["raw_expression_missing_fraction_by_gene"][gene] <= 0.05
    ]
    if len(common_genes) < 45:
        raise ValueError(
            f"Only {len(common_genes)} PAM50 genes passed the cross-cohort "
            "95% coverage rule; at least 45 are required."
        )
    gene_indices = [molecular_model.PAM50_GENES.index(gene) for gene in common_genes]
    tcga_raw = tcga["expression_values"][:, gene_indices]
    metabric_raw = metabric["expression_values"][:, gene_indices]

    tcga_rank = molecular_model.rank_normalise_expression(
        tcga_raw,
        common_genes,
        minimum_observed_fraction=0.95,
    )
    metabric_rank = molecular_model.rank_normalise_expression(
        metabric_raw,
        common_genes,
        minimum_observed_fraction=0.95,
    )
    tcga_scores, metabric_scores, pca_metadata = molecular_model.fit_tcga_pca_projection(
        tcga_rank.values,
        metabric_rank.values,
        common_genes,
        n_components=n_components,
        random_state=seed,
    )
    results = _fit_cox_model(
        tcga["patient_rows"],
        metabric["patient_rows"],
        tcga_scores,
        metabric_scores,
        penalizer=cox_penalizer,
        n_bootstrap=n_bootstrap,
        seed=seed,
    )
    subtype_transfer = molecular_model.evaluate_intrinsic_subtype_transfer(
        tcga_rank.values,
        tcga["patient_ids"],
        tcga["patient_subtype_labels"],
        metabric_rank.values,
        metabric["patient_ids"],
        metabric["patient_subtype_labels"],
        c_value=0.1,
        n_bootstrap=n_bootstrap,
        seed=seed,
    )

    cohort_metadata = []
    for cohort, rank_report in (
        (tcga, tcga_rank),
        (metabric, metabric_rank),
    ):
        cohort_metadata.append(
            {
                "label": cohort["label"],
                "study_id": cohort["study_id"],
                "profile_id": cohort["profile_id"],
                "profile_name": cohort["profile_name"],
                "profile_description": cohort["profile_description"],
                "expression_samples_in_profile_list": cohort["n_expression_samples_in_list"],
                "primary_tumour_samples": cohort["n_primary_tumour_samples"],
                "expression_samples_used": cohort["n_expression_samples_used"],
                "patients_in_survival_expression_analysis": cohort[
                    "n_unique_primary_tumour_patients_with_os"
                ],
                "death_events": cohort["n_os_events"],
                "censored": cohort["n_censored"],
                "sample_type_filter": cohort["sample_type_filter"],
                "sample_type_counts": cohort["sample_type_counts"],
                "duplicate_primary_samples_excluded": cohort[
                    "duplicate_primary_samples_excluded"
                ],
                "samples_without_patient_link": cohort["samples_without_patient_link"],
                "expression_fetch_qc": cohort["expression_fetch_qc"],
                "samples_with_at_least_one_imputed_gene": rank_report.samples_with_imputation,
                "genes": len(common_genes),
                "missing_fraction_by_gene": {
                    gene: rank_report.missing_fraction_by_gene[gene]
                    for gene in common_genes
                },
                "subtype_annotation_audit": cohort["subtype_annotation_audit"],
            }
        )

    return {
        "analysis": "PAM50_panel_rank_normalised_TCGA_to_METABRIC_survival_validation",
        "status": "completed",
        "date": "2026-11-11",
        "panel": {
            "name": "PAM50 50-gene intrinsic-subtype research panel",
            "panel_size": len(molecular_model.PAM50_GENES),
            "common_genes_used": common_genes,
            "common_gene_count": len(common_genes),
            "resolved_hugo_symbols": resolved_symbols,
            "gene_identifier_source": (
                "cBioPortal gene endpoint; aliases resolved to current HUGO symbols"
            ),
        },
        "cohorts": cohort_metadata,
        "harmonisation": {
            "method": (
                "Per-gene within-cohort mid-rank inverse-normal transformation. "
                "Applied independently to TCGA and METABRIC expression values."
            ),
            "uses_external_clinical_outcomes": False,
            "uses_external_subtype_labels": False,
            "external_expression_distribution_used": True,
            "pca_fit_only_on_tcga": True,
            "pca": pca_metadata,
            "transductive": True,
            "individual_sample_deployment_supported": False,
            "limitations": [
                (
                    "Cohort-wise rank normalisation uses the full unlabeled expression "
                    "distribution of each cohort."
                ),
                (
                    "Gene-level cBioPortal matrices are used; the underlying METABRIC "
                    "probe-level annotation and duplicate-probe resolution are not "
                    "rederived here."
                ),
                (
                    "Rank normalisation reduces scale differences but cannot remove "
                    "all platform-specific gene measurement bias."
                ),
                (
                    "PCA captures variance in the TCGA feature matrix and is not a "
                    "validated intrinsic-subtype classifier."
                ),
            ],
        },
        "survival_endpoint": {
            "name": "overall survival",
            "fields": ["OS_MONTHS", "OS_STATUS"],
            "event_encoding": "1:DECEASED=1; 0:LIVING=0",
            "time_unit": "months",
            "time_origin": (
                "cBioPortal documents OS_MONTHS as months from a study-defined starting point; "
                "source-specific endpoint provenance remains under review."
            ),
            "nonpositive_followup": "excluded; never replaced with arbitrary positive time",
        },
        "model_results": results,
        "subtype_transfer_results": subtype_transfer,
        "prespecification": {
            "secondary_subtype_task": (
                "TCGA-trained, METABRIC-tested five-class logistic classifier "
                "on rank-normalised fixed 50-gene features"
            ),
            "subtype_classifier_C": 0.1,
            "subtype_shared_classes": list(molecular_model.INTRINSIC_SUBTYPE_CLASSES),
            "primary_model": "age + first five TCGA-fitted PAM50-panel principal components",
            "age_only_comparator": "age only, fitted on the same complete-case TCGA patients",
            "pca_components": int(n_components),
            "cox_penalizer": float(cox_penalizer),
            "bootstrap_replicates": int(n_bootstrap),
            "bootstrap_seed": int(seed),
            "external_outcomes_used_for_feature_construction_or_tuning": False,
            "external_subtype_labels_used_for_feature_construction_or_tuning": False,
        },
        "data_governance": {
            "patient_level_rows_written": False,
            "sample_ids_written": False,
            "patient_ids_written": False,
            "raw_expression_values_written": False,
            "aggregate_metrics_only": True,
        },
        "interpretation": (
            "Exploratory molecular association and external transport evaluation. "
            "It is not a commercial PAM50/Prosigna classifier, clinical prognostic test, "
            "treatment-response model, or validation of the spatial/mechanistic digital twin."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default="results/molecular_validation/pam50_tcga_metabric_survival.json",
    )
    parser.add_argument("--components", type=int, default=N_COMPONENTS)
    parser.add_argument("--bootstrap", type=int, default=250)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--cox-penalizer", type=float, default=COX_PENALIZER)
    args = parser.parse_args()
    if args.bootstrap < 0:
        parser.error("--bootstrap must be non-negative")
    report = run_pam50_molecular_validation(
        n_components=args.components,
        n_bootstrap=args.bootstrap,
        seed=args.seed,
        cox_penalizer=args.cox_penalizer,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for cohort in report["cohorts"]:
        print(
            cohort["label"],
            "patients:", cohort["patients_in_survival_expression_analysis"],
            "events:", cohort["death_events"],
            "genes:", cohort["genes"],
        )
    result = report["model_results"]
    print("METABRIC age-only C-index:", result["age_only_reference"]["harrell_c_index"])
    print("METABRIC molecular C-index:", result["molecular_model"]["harrell_c_index"])
    print("Delta C-index:", result["molecular_model"]["delta_c_index_vs_age_only"])
    print(
        "Subtype transfer balanced accuracy:",
        report["subtype_transfer_results"]["external_evaluation"]["balanced_accuracy"],
    )
    print(
        "Subtype transfer macro-F1:",
        report["subtype_transfer_results"]["external_evaluation"]["macro_f1"],
    )
    print("Aggregate molecular report:", output)


if __name__ == "__main__":
    main()
