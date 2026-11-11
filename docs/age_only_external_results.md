# Age-only external survival validation — results

**Result:** completed successfully  
**Date:** 11 November 2026  
**Workflow:** [GitHub Actions run 42980553346](https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/42980553346)  
**Evaluated source commit:** f115d52f9d25e188101a5b071cebe67c84938e07  
**Aggregate JSON:** results/survival_validation/age_only_external_metrics_summary.json

## Cohort flow

| Cohort | Patients in analysis | Death events | Censored |
|---|---:|---:|---:|
| TCGA-BRCA development | 1,069 | 151 | 918 |
| METABRIC external evaluation | 1,979 | 1,143 | 836 |

The development cohort began with 1,082 RNA-seq expression sample-list members; 13 non-positive OS times were excluded. The external cohort began with 1,980 expression sample-list members; one non-positive OS time was excluded. No time was replaced with an arbitrary positive value.

## Model

A Cox proportional-hazards model was fitted in TCGA using age at diagnosis only. Age was standardized using the TCGA mean (58.27 years) and standard deviation (13.16 years), and those values were reused unchanged in METABRIC. Penalizer = 0.1. Stage was deliberately excluded because cross-cohort stage semantics remain unresolved. Molecular expression was not used.

The age hazard ratio per one TCGA standard deviation was **1.277** (95% CI 1.129–1.445; p = 9.94e-5). The training-cohort rank-based proportional-hazards test p-value was 0.284. A non-significant test does not prove the proportional-hazards assumption.

## External performance

### Global discrimination

| Metric | METABRIC result |
|---|---:|
| Harrell C-index | **0.5882** |
| Patient bootstrap 95% CI | 0.5683–0.6065 |
| Bootstrap replicates | 250 |

This is weak discrimination. It is a reference baseline, not a strong prognostic model.

### Time-dependent discrimination and prediction error

All time-dependent metrics use inverse-probability-of-censoring weighting (IPCW). The censoring Kaplan-Meier estimator is calculated within METABRIC for metric correction only; model fitting and individual predictions remain based on TCGA.

| Horizon | Cases | Controls | IPCW AUC (95% CI) | IPCW Brier (95% CI) |
|---|---:|---:|---:|---:|
| 12 months | 30 | 1,934 | 0.5783 (0.4602–0.6844) | 0.0150 (0.0093–0.0204) |
| 36 months | 223 | 1,717 | 0.5151 (0.4740–0.5613) | 0.1014 (0.0904–0.1119) |
| 60 months | 427 | 1,489 | 0.5395 (0.5037–0.5719) | 0.1713 (0.1614–0.1818) |

The one-year and three-year AUC intervals include 0.5. The five-year interval is only slightly above 0.5. Overall, age alone provides little separation of short- and medium-term outcomes in this external cohort.

The sparse-grid integrated Brier score over time zero and 12, 36 and 60 months was **0.07935**. This is a trapezoidal approximation on the prespecified sparse grid, not a dense-grid integrated Brier score.

### Calibration summary

The table compares overall mean predicted survival with the IPCW estimate of observed survival. The gap is observed minus predicted.

| Horizon | Mean predicted survival | IPCW observed survival | Gap (95% bootstrap CI) |
|---|---:|---:|---:|
| 12 months | 0.9783 | 0.9847 | +0.0064 (+0.0009 to +0.0124) |
| 36 months | 0.8916 | 0.8858 | -0.0058 (-0.0192 to +0.0090) |
| 60 months | 0.8019 | 0.7801 | -0.0218 (-0.0398 to -0.0052) |

At 60 months, the model overestimates survival by about 2.2 percentage points on average in this external cohort. The decile calibration data in the machine-readable artifact also show heterogeneity across predicted-survival groups; a single overall gap does not capture all calibration errors.

## Interpretation

The full workflow passed repository tests, lint checks and the external cohort evaluation. It provides a reproducible age-only reference model, but the external discrimination is weak. This result argues against treating age-only prediction as a useful standalone model and establishes a baseline that future molecular and mechanistic models must improve on under the same locked evaluation design.

This run does **not** validate the ONCO-DYNAMICS-X digital twin. It uses no expression, spatial imaging, cell-graph, PDE, agent-based, PINN, neural-operator or treatment-response features. It does not establish clinical utility, patient-level prognosis, treatment benefit or a causal effect of age.

## Method limitations

- The IPCW correction uses a marginal censoring model and assumes independent censoring within METABRIC; informative or covariate-dependent censoring is not modelled.
- The stage predictor was intentionally excluded due to the unresolved difference between TCGA patient-level AJCC pathologic stage and METABRIC sample-level generic tumour stage.
- The C-index and time-dependent AUC summarize discrimination, not calibration or net clinical benefit.
- The Brier score is sensitive to prediction calibration and the stated censoring assumptions.
- Integrated Brier score is based on only three follow-up horizons plus time zero.
- This is a retrospective public-cohort transport experiment, not prospective or clinical validation.

## Next decision

Keep this age-only model as the primary clinical reference comparator. Do not reintroduce stage into the primary cross-cohort model until the original METABRIC stage definition, assessment timepoint and category 0 meaning are verified. The next molecular workstream must lock a fixed gene/signature panel and cross-platform mapping before any METABRIC outcome-based evaluation, and must be compared against this baseline with calibration as well as discrimination.
