# Independent cohort and clinical-outcome validation protocol

## Purpose

Move ONCO-DYNAMICS X beyond the initial TCGA-BRCA tumour-versus-normal pipeline without treating specimen classification as evidence of prognosis, treatment response, or a validated digital twin.

This document is a prespecified analysis protocol. It does not claim that an independent cohort has already been downloaded, harmonized, or analysed.

## 1. Separate the research questions

Keep these endpoints distinct:

1. **Specimen-label engineering benchmark:** distinguish primary tumour from solid-tissue normal. This checks the expression-data pipeline only.
2. **Prognostic analysis:** estimate association between a molecular score and a defined time-to-event outcome in patients with breast cancer.
3. **Treatment-response analysis:** estimate response conditional on a documented treatment and an explicit response definition. A prognostic association is not evidence of treatment benefit.
4. **Spatial/mechanistic validation:** compare spatial model fields or cell-interaction claims with spatially resolved tissue measurements. Bulk RNA-seq cannot validate cell coordinates, cell contacts, oxygen gradients, or drug penetration.

Do not pool these questions into one model or report one as evidence for another.

## 2. Candidate independent cohort: METABRIC

METABRIC is a candidate external breast-cancer cohort, accessible through cBioPortal. Before analysis, perform and archive a data audit covering:

- Study release/version and access terms.
- Number of unique patients and samples; duplicate and multiple-sample handling.
- Molecular assay/platform and expression units.
- Clinical variable dictionary and missingness.
- Exact definitions and availability of overall survival, disease-specific survival, relapse-free interval or other event endpoints.
- Time origin, event indicator, censoring convention, follow-up units and implausible values.
- Treatment exposure fields and response labels, if any.
- Patient/sample identifiers linking molecular and clinical records.

The candidate cohort is accepted only if the endpoint and molecular data are linked at the patient level and the access terms permit the intended use. Record an exclusion reason if any requirement fails. Do not assume that all candidate endpoints or treatment variables are available in a given study release.

## 3. Freeze the TCGA development analysis

Before examining external-cohort outcome performance:

- Freeze the TCGA development cohort, inclusion rules, preprocessing, feature list, model family and hyperparameters.
- Preserve patient-level grouping; no patient may contribute samples to more than one partition.
- Record the exact Git commit, GDC file UUIDs, MD5 checksums, query filters, software versions and random seeds.
- Keep all feature selection, scaling, imputation and tuning inside the development resampling process.
- Do not select features or thresholds using METABRIC outcomes.
- Preserve a simple clinical-only or null baseline where appropriate.

The current paired tumour/normal task is not the development target for survival prediction. A separate, endpoint-specific model must be specified and trained using eligible tumour patients and outcome labels.

## 4. Molecular harmonization

TCGA STAR-count/TPM measurements and METABRIC microarray measurements are not interchangeable.

- Confirm platform, probe annotation, gene identifier version and expression scale.
- Map probes to genes using a documented annotation release; handle many-to-one mappings explicitly.
- Restrict to a defensible common feature space before scoring.
- Apply transformations and normalization using a protocol that does not use external outcome labels.
- Do not concatenate cohorts and randomly split samples; that can leak cohort/platform signatures.
- Prefer a frozen TCGA-derived model evaluated directly in METABRIC after justified platform harmonization. If direct transport is not technically defensible, report that external validation is not yet possible rather than silently retraining on the test cohort.
- Include platform/batch diagnostics and a simple cohort-origin classifier as a warning check, not as proof of biological validity.

## 5. Prognostic endpoint analysis

For each candidate endpoint, write a separate endpoint specification before fitting:

- Event definition and source field.
- Time origin and time unit.
- Event indicator and censoring rules.
- Eligibility and exclusions.
- Minimum follow-up and missingness handling.

Report the number of patients, events, censored observations and median follow-up. Use Kaplan-Meier curves as descriptive summaries and Cox proportional-hazards models only after checking assumptions. Report hazard ratios with confidence intervals, calibration at prespecified time horizons where feasible, and discrimination such as Harrell's C-index or time-dependent AUC with appropriate uncertainty.

Check proportional hazards, influential observations, missing-data sensitivity and subgroup stability. Avoid dichotomizing continuous risk scores unless a threshold was fixed in the development cohort. If the event count is inadequate, report descriptive results and do not fit an unstable multivariable model.

A statistically significant association alone is not clinical utility. Assess calibration and incremental value over available clinical covariates before making stronger claims.

## 6. Treatment-response analysis

Proceed only if the cohort has patient-linked treatment exposure and a valid response endpoint with a documented definition and assessment time.

- Specify treatment, line/context, response criteria, assessment window and missing/unevaluable handling.
- Distinguish pathological complete response, objective response, recurrence and survival; they are not interchangeable.
- Account for confounding by indication and treatment selection. Observational treatment comparisons do not establish causal treatment benefit.
- Do not infer checkpoint-inhibitor response from general survival or from tumour-versus-normal labels.
- If treatment/response linkage is absent or too sparse, mark this analysis as unavailable for that cohort.

## 7. Validation hierarchy and reporting

Use the following labels consistently:

- **Software verified:** tests, deterministic runs, input validation and reproducibility checks pass.
- **Numerically verified:** solver convergence/stability and independent numerical comparisons pass.
- **Internally evaluated:** patient-grouped development evaluation is complete.
- **Externally evaluated:** a frozen model is tested on a separately sourced, compatible cohort.
- **Biologically supported:** the relevant model component agrees with modality-appropriate biological measurements.
- **Clinically validated:** requires substantially more evidence and appropriate prospective/clinical validation; this project must not use this label based on retrospective public-cohort analysis alone.

Report failed analyses, missingness, cohort flow, exclusions, uncertainty and negative results. Do not report unrun analyses as completed.

## 8. Data governance

Use only data whose current access terms permit the analysis. Keep controlled, restricted or identifiable participant-level data out of this public repository. Commit code, schemas, aggregate metrics and non-sensitive manifests only where permitted. Do not commit raw expression matrices or clinical tables. Record data versions and checksums so an authorized user can reproduce the work locally.

## Candidate source

- METABRIC study page in cBioPortal: https://www.cbioportal.org/study/summary?id=brca_metabric
- cBioPortal documentation: https://docs.cbioportal.org/

The study page is a candidate discovery point, not a guarantee that every molecular field, endpoint or treatment variable needed by this protocol is available.
