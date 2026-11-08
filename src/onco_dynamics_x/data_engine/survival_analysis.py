"""Exploratory cross-cohort clinical survival baseline.

This module fits a simple age + stage Cox model in TCGA-BRCA and evaluates
the frozen model in METABRIC. It is not a clinical prediction tool.
"""
from __future__ import annotations

import json
import re
from collections import Counter
from typing import Any

from onco_dynamics_x.data_engine.cbioportal import CBioPortalClient


def normalise_stage(value: Any) -> int | None:
    """Map common stage labels to broad anatomic groups I-IV."""
    if value is None:
        return None
    text = re.sub(r"[^A-Z0-9]", "", str(value).upper())
    if not text or text in {"NA", "NAN", "UNKNOWN", "NOTREPORTED", "NOTAVAILABLE"}:
        return None
    text = text.replace("STAGE", "")
    if text.startswith("IV") or text in {"4", "4A", "4B", "4C"}:
        return 4
    if text.startswith("III") or text in {"3", "3A", "3B", "3C"}:
        return 3
    if text.startswith("II") or text in {"2", "2A", "2B", "2C"}:
        return 2
    if text.startswith("I") or text in {"1", "1A", "1B", "1C"}:
        return 1
    return None


def parse_os_event(value: Any) -> int | None:
    """cBioPortal convention: 1 = death/event; 0 = living/censored."""
    text = str(value or "").strip().upper()
    if text.startswith("1:"):
        return 1
    if text.startswith("0:"):
        return 0
    return None


def _attribute_values(records, attribute_id, entity_key):
    values = {}
    for row in records:
        if row.get("clinicalAttributeId") != attribute_id:
            continue
        entity_id = row.get(entity_key)
        value = row.get("value")
        if entity_id and value is not None and str(value).strip() not in {"", "NA", "NaN"}:
            values[entity_id] = value
    return values


def build_survival_records(
    *,
    patients: list[dict[str, Any]],
    samples: list[dict[str, Any]],
    expression_sample_ids: list[str],
    patient_clinical: list[dict[str, Any]],
    sample_clinical: list[dict[str, Any]],
    age_attribute: str,
    stage_attribute: str,
    stage_level: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Assemble a de-identified-in-memory patient table and aggregate QC."""
    patient_ids = {row.get("patientId") for row in patients if row.get("patientId")}
    sample_to_patient = {
        row.get("sampleId"): row.get("patientId")
        for row in samples
        if row.get("sampleId") and row.get("patientId") in patient_ids
    }
    age_values = _attribute_values(patient_clinical, age_attribute, "patientId")
    time_values = _attribute_values(patient_clinical, "OS_MONTHS", "patientId")
    status_values = _attribute_values(patient_clinical, "OS_STATUS", "patientId")
    if stage_level == "PATIENT":
        stage_values = _attribute_values(patient_clinical, stage_attribute, "patientId")
        stage_key = "patientId"
    elif stage_level == "SAMPLE":
        stage_values = _attribute_values(sample_clinical, stage_attribute, "sampleId")
        stage_key = "sampleId"
    else:
        raise ValueError("stage_level must be PATIENT or SAMPLE")

    rows = []
    seen_patients = set()
    counts = Counter()
    stage_counts = Counter()
    event_counts = Counter()
    for sample_id in expression_sample_ids:
        patient_id = sample_to_patient.get(sample_id)
        if not patient_id:
            counts["expression_samples_without_patient_link"] += 1
            continue
        if patient_id in seen_patients:
            counts["duplicate_expression_samples_after_patient_grouping"] += 1
            continue
        seen_patients.add(patient_id)
        counts["expression_samples_linked"] += 1

        if patient_id not in time_values or patient_id not in status_values:
            counts["missing_os_pair"] += 1
            continue
        try:
            duration = float(time_values[patient_id])
        except (TypeError, ValueError):
            counts["unparseable_os_time"] += 1
            continue
        if duration <= 0:
            counts["nonpositive_os_time"] += 1
            continue
        event = parse_os_event(status_values[patient_id])
        if event is None:
            counts["unrecognized_os_status"] += 1
            continue
        try:
            age = float(age_values[patient_id])
        except (KeyError, TypeError, ValueError):
            counts["missing_or_invalid_age"] += 1
            continue
        stage_value = (
            stage_values.get(patient_id)
            if stage_key == "patientId"
            else stage_values.get(sample_id)
        )
        stage = normalise_stage(stage_value)
        if stage is None:
            counts["missing_or_unrecognized_stage"] += 1
            continue
        rows.append(
            {
                "patient_id": patient_id,
                "duration_months": duration,
                "event": event,
                "age_years": age,
                "stage_group": stage,
            }
        )
        stage_counts[str(stage)] += 1
        event_counts["deceased" if event else "censored"] += 1

    ages = [row["age_years"] for row in rows]
    durations = [row["duration_months"] for row in rows]
    summary = {
        "expression_sample_list_count": len(expression_sample_ids),
        "unique_expression_patients_seen": len(seen_patients),
        "patients_retained_for_model": len(rows),
        "event_counts": dict(sorted(event_counts.items())),
        "stage_group_counts": dict(sorted(stage_counts.items())),
        "exclusion_counts": dict(sorted(counts.items())),
        "age_range_years": (
            [round(min(ages), 3), round(max(ages), 3)]
            if ages else None
        ),
        "followup_range_months": (
            [round(min(durations), 4), round(max(durations), 4)]
            if durations else None
        ),
    }
    return rows, summary


def _design_matrix(frame, age_mean: float, age_sd: float):
    import pandas as pd

    design = pd.DataFrame(index=frame.index)
    design["age_z"] = (frame["age_years"] - age_mean) / age_sd
    for stage in (2, 3, 4):
        design[f"stage_{stage}"] = (frame["stage_group"] == stage).astype(float)
    return design


def run_cross_cohort_cox(
    tcga_records: list[dict[str, Any]],
    metabric_records: list[dict[str, Any]],
    penalizer: float = 0.1,
) -> dict[str, Any]:
    """Fit age+stage Cox in TCGA and evaluate once on untouched METABRIC."""
    import numpy as np
    import pandas as pd

    tcga = pd.DataFrame(tcga_records)
    metabric = pd.DataFrame(metabric_records)
    try:
        from lifelines import CoxPHFitter
        from lifelines.statistics import proportional_hazard_test
        from lifelines.utils import concordance_index
    except ImportError as exc:
        raise RuntimeError(
            "Install the optional survival dependencies with pip install -e '.[survival]'."
        ) from exc

    if len(tcga) < 100 or int(tcga["event"].sum()) < 30:
        raise ValueError("TCGA development cohort is too small for this baseline.")
    if len(metabric) < 100 or int(metabric["event"].sum()) < 30:
        raise ValueError("METABRIC external cohort is too small for this baseline.")

    age_mean = float(tcga["age_years"].mean())
    age_sd = float(tcga["age_years"].std(ddof=0))
    if not np.isfinite(age_sd) or age_sd == 0:
        raise ValueError("TCGA age variance is zero or invalid.")

    train_design = _design_matrix(tcga, age_mean, age_sd)
    test_design = _design_matrix(metabric, age_mean, age_sd)
    features = [
        column for column in train_design.columns
        if train_design[column].nunique() > 1
    ]
    if not features:
        raise ValueError("No variable predictors remain after design-matrix checks.")

    train = train_design[features].copy()
    train["duration_months"] = tcga["duration_months"].to_numpy()
    train["event"] = tcga["event"].to_numpy()
    test = test_design[features].copy()
    test["duration_months"] = metabric["duration_months"].to_numpy()
    test["event"] = metabric["event"].to_numpy()

    model = CoxPHFitter(penalizer=penalizer)
    model.fit(train, duration_col="duration_months", event_col="event")
    test_hazard = model.predict_partial_hazard(test[features]).to_numpy().reshape(-1)
    test_durations = test["duration_months"].to_numpy()
    test_events = test["event"].to_numpy()
    external_c_index = float(
        concordance_index(test_durations, -test_hazard, test_events)
    )

    age_train = train[["age_z", "duration_months", "event"]].copy()
    age_test = test[["age_z"]].copy()
    age_model = CoxPHFitter(penalizer=penalizer)
    age_model.fit(age_train, duration_col="duration_months", event_col="event")
    age_hazard = age_model.predict_partial_hazard(age_test).to_numpy().reshape(-1)
    age_only_c_index = float(
        concordance_index(test_durations, -age_hazard, test_events)
    )

    def bootstrap_ci(predicted_hazard):
        rng = np.random.default_rng(42)
        estimates = []
        n = len(test_durations)
        for _ in range(500):
            indices = rng.integers(0, n, size=n)
            try:
                estimate = concordance_index(
                    test_durations[indices],
                    -predicted_hazard[indices],
                    test_events[indices],
                )
            except (ZeroDivisionError, ValueError):
                continue
            if np.isfinite(estimate):
                estimates.append(float(estimate))
        if not estimates:
            return None
        return [
            float(np.quantile(estimates, 0.025)),
            float(np.quantile(estimates, 0.975)),
        ]

    ph_test = proportional_hazard_test(
        model, train, time_transform="rank"
    ).summary
    coefficients = []
    for variable, row in model.summary.iterrows():
        coefficients.append(
            {
                "variable": str(variable),
                "coefficient": float(row["coef"]),
                "hazard_ratio": float(row["exp(coef)"]),
                "ci_lower_hazard_ratio": float(row["exp(coef) lower 95%"]),
                "ci_upper_hazard_ratio": float(row["exp(coef) upper 95%"]),
                "p_value": float(row["p"]),
            }
        )
    return {
        "model": "Cox proportional hazards; age + broad stage groups",
        "training_cohort": {
            "patients": len(tcga),
            "events": int(tcga["event"].sum()),
            "censored": int((1 - tcga["event"]).sum()),
        },
        "external_cohort": {
            "patients": len(metabric),
            "events": int(metabric["event"].sum()),
            "censored": int((1 - metabric["event"]).sum()),
            "harrell_c_index": external_c_index,
            "harrell_c_index_bootstrap_95ci": bootstrap_ci(test_hazard),
            "age_only_harrell_c_index": age_only_c_index,
            "age_only_harrell_c_index_bootstrap_95ci": bootstrap_ci(age_hazard),
        },
        "configuration": {
            "age_mean_from_tcga": age_mean,
            "age_sd_from_tcga": age_sd,
            "penalizer": penalizer,
            "bootstrap_replicates": 500,
            "bootstrap_seed": 42,
            "features": features,
            "stage_reference_group": 1,
            "stage_mapping": "Broad AJCC stage I-IV; unknown/0 excluded",
        },
        "coefficients": coefficients,
        "proportional_hazards_test_p_values": {
            str(index): float(row["p"])
            for index, row in ph_test.iterrows()
        },
        "interpretation": {
            "external_evaluation_performed": True,
            "molecular_model": False,
            "clinical_validation": False,
            "treatment_response_prediction": False,
            "note": (
                "Exploratory clinical-only transport test. METABRIC outcomes were "
                "not used for fitting or tuning. Stage harmonization and proportional "
                "hazards assumptions require review before interpretation."
            ),
        },
    }


def write_survival_report(report: dict[str, Any], path: str) -> None:
    from pathlib import Path

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def fetch_survival_cohort(
    client: CBioPortalClient,
    study_id: str,
    sample_list_id: str,
    age_attribute: str,
    stage_attribute: str,
    stage_level: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Fetch one public cBioPortal cohort and assemble its endpoint/covariates."""
    patients = client.get_all(
        f"studies/{study_id}/patients", {"projection": "SUMMARY"}
    )
    samples = client.get_all(
        f"studies/{study_id}/samples", {"projection": "SUMMARY"}
    )
    patient_clinical = client.get_all(
        f"studies/{study_id}/clinical-data", {"clinicalDataType": "PATIENT"}
    )
    sample_clinical = client.get_all(
        f"studies/{study_id}/clinical-data", {"clinicalDataType": "SAMPLE"}
    )
    sample_list = client.get_json(f"sample-lists/{sample_list_id}")
    return build_survival_records(
        patients=patients,
        samples=samples,
        expression_sample_ids=sample_list.get("sampleIds", []),
        patient_clinical=patient_clinical,
        sample_clinical=sample_clinical,
        age_attribute=age_attribute,
        stage_attribute=stage_attribute,
        stage_level=stage_level,
    )
