# TCGA-BRCA paired tumour/normal baseline: expanded 48-patient run

## Run provenance

- Repository: https://github.com/culprit-k/ONCO-DYNAMICS-X
- Workflow run: https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/42452284524
- Workflow: `TCGA-BRCA paired expression baseline`
- Workflow run ID: `42452284524`
- Run commit: `2e5a4bd2bdc7a06c7b5a7142c6d151326416c269`
- Conclusion: **success**
- Artifact: `tcga-brca-paired-baseline-report` (artifact ID `13605679881`)
- Artifact contains selected metadata/manifests and metrics JSON, not raw expression files or the expression matrix.

## Cohort and pipeline

The completed run selected 48 patients with one primary-tumour and one solid-tissue-normal expression sample per patient:

| Item | Count |
|---|---:|
| Patients | 48 |
| Expression samples | 96 |
| Primary tumour samples | 48 |
| Solid tissue normal samples | 48 |
| Holdout training patients | 36 |
| Holdout test patients | 12 |
| Holdout training samples | 72 |
| Holdout test samples | 24 |

Source: NCI GDC, project `TCGA-BRCA`, `STAR - Counts` workflow, `tpm_unstranded` column, open-access files. The workflow completed file downloads with checksum verification, assembled the expression matrix, ran the patient-grouped evaluation and uploaded its report artifact.

The model pipeline used `log1p(TPM)`, zero-variance filtering, ANOVA top-100 feature selection, standardization and class-balanced logistic regression. Feature selection and scaling were fitted inside the training folds. The holdout and five-fold cross-validation both grouped by patient.

## Metrics as reported by the workflow

| Metric | Result |
|---|---:|
| Holdout balanced accuracy | 1.000 |
| Holdout ROC-AUC | 1.000 |
| Five-fold grouped CV balanced accuracy (mean ± SD) | 1.000 ± 0.000 |
| Five-fold grouped CV ROC-AUC (mean ± SD) | 1.000 ± 0.000 |

Holdout confusion matrix (normal first, tumour second):

| Actual / predicted | Normal | Tumour |
|---|---:|---:|
| Normal | 12 | 0 |
| Tumour | 0 | 12 |

These are the observed metrics from this workflow run, not estimates of real-world diagnostic performance.

## Interpretation and important limitations

This is an **engineering benchmark for distinguishing known specimen labels** in a small paired cohort. It does not predict future cancer, recurrence, survival, metastasis, treatment response or benefit. It does not validate the PDE, PINN, agent-based, graph, immune, hypoxia or drug-diffusion components of ONCO-DYNAMICS X.

Perfect scores deserve additional scrutiny rather than celebration as clinical performance. The two classes are biologically distinct tissue specimen types, and the task is much easier than patient-level clinical decision-making. The cohort is small, selected for paired tumour/normal availability, and comes from one project and assay workflow. Patient grouping prevents a patient's paired specimens from crossing the train/test boundary, but does not remove all possible selection, technical or cohort-specific effects.

Before using any molecular feature as a biomarker, perform independent cohort evaluation, feature stability analysis, platform/batch diagnostics and appropriate endpoint-specific validation. Do not interpret selected logistic-regression coefficients as causal genes or therapeutic targets.

## Next validation gates

1. Audit METABRIC/cBioPortal molecular and clinical records, access terms, endpoint definitions and patient/sample linkage.
2. Define a separate tumour-patient prognostic question with time-to-event and censoring data; the current tumour/normal classifier is not a survival model.
3. Freeze any development model before external evaluation and assess platform harmonization between TCGA RNA-seq and METABRIC microarray data.
4. Report endpoint-specific discrimination, calibration, uncertainty and clinical-covariate baselines only if event counts and data quality support them.
5. Keep treatment response separate; analyse it only if treatment exposure and a valid response endpoint are actually linked.
6. Validate spatial and mechanistic model components with spatially resolved data, not bulk expression alone.

See [the independent validation protocol](independent_validation_protocol.md) for the prespecified cohort audit and analysis requirements.
