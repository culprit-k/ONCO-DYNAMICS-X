# Exploratory TCGA-to-METABRIC clinical survival baseline

## Provenance

- Repository: https://github.com/culprit-k/ONCO-DYNAMICS-X
- Workflow run: https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/42454365253
- Workflow conclusion: success
- Run commit: 7df0c3913cadf9d132680b0a0730982d6b6cfc95
- Artifact: tcga-metabric-clinical-survival-baseline
- Model family: Cox proportional hazards
- Endpoint: overall survival (OS_MONTHS, OS_STATUS)
- Development cohort: TCGA-BRCA PanCancer Atlas
- External evaluation cohort: METABRIC

The workflow retrieves public clinical metadata in memory and uploads only aggregate QC and model metrics. It does not upload patient IDs, sample IDs, row-level clinical data or expression matrices.

## Analysis design

A deliberately simple clinical-only model was fitted in TCGA and evaluated once in METABRIC.

- Predictors: age at diagnosis and broad stage group I-IV.
- TCGA age field: AGE.
- TCGA stage field: AJCC_PATHOLOGIC_TUMOR_STAGE.
- METABRIC age field: AGE_AT_DIAGNOSIS.
- METABRIC stage field: sample-level TUMOR_STAGE, linked to the patient through the expression sample.
- Endpoint: OS_MONTHS and OS_STATUS.
- Event convention: 1 = deceased; 0 = living/censored, following the cBioPortal OS field convention.
- Non-positive follow-up times were excluded and counted explicitly; they were not silently changed.
- Age was standardized using the TCGA development cohort mean and standard deviation only.
- Stage I was the reference group; stages II-IV were one-hot encoded.
- Cox penalizer: 0.1, fixed before external evaluation.
- No METABRIC outcome was used for feature selection, fitting or tuning.
- External C-index uncertainty was estimated with 500 patient-level bootstrap resamples (seed 42).

This is an exploratory transport test of a small clinical baseline. It is not a molecular model and does not validate the ONCO-DYNAMICS X spatial digital twin.

## Cohort flow

| Item | TCGA development | METABRIC external |
|---|---:|---:|
| Expression samples queried | 1,082 | 1,980 |
| Expression samples linked to a patient | 1,082 | 1,980 |
| Patients retained for Cox model | 1,050 | 1,454 |
| Death events | 143 | 821 |
| Living/censored | 907 | 633 |
| Non-positive OS time excluded | 13 | 1 |
| Missing/unrecognized stage excluded | 19 | 525 |

METABRIC has substantial stage missingness in the expression subset: 525 of 1,980 patients were excluded for missing or unrecognized stage. The final external cohort is therefore a complete-case subset, not the entire METABRIC expression cohort.

Stage distribution in the retained data:

| Broad stage | TCGA | METABRIC |
|---|---:|---:|
| I | 180 | 501 |
| II | 606 | 825 |
| III | 245 | 118 |
| IV | 19 | 10 |

Stage IV is sparse in both cohorts, especially in the external cohort (10 patients). Its estimated hazard ratio is consequently imprecise.

## External discrimination

| Model | METABRIC Harrell C-index | Bootstrap 95% interval |
|---|---:|---:|
| Age only | 0.5914 | 0.5700–0.6174 |
| Age + broad stage | 0.6161 | 0.5923–0.6411 |

The age+stage model's C-index is higher than the age-only comparator in this evaluation, but the intervals overlap. Do not claim that stage provides a statistically established improvement based on this comparison alone. A C-index around 0.62 is modest discrimination, not sufficient evidence of clinical utility.

## TCGA-fitted coefficients

Hazard ratios are relative to stage I, with age expressed per one TCGA standard deviation (13.16 years).

| Predictor | Hazard ratio | 95% CI | p-value |
|---|---:|---:|---:|
| Age (per 1 SD) | 1.273 | 1.121–1.444 | 0.000188 |
| Stage II vs I | 0.983 | 0.747–1.293 | 0.902 |
| Stage III vs I | 1.589 | 1.163–2.170 | 0.00365 |
| Stage IV vs I | 5.854 | 3.283–10.440 | 2.12e-9 |

These are observational associations in the TCGA development data, not causal effects.

## Proportional-hazards diagnostics

Rank-based proportional-hazards test p-values:

| Predictor | p-value |
|---|---:|
| Age | 0.333 |
| Stage II | 0.422 |
| Stage III | 0.0557 |
| Stage IV | 0.0744 |

No test crossed the conventional 0.05 threshold, but stage III and IV are close to it. This does not prove the proportional-hazards assumption. The sparse stage IV group and borderline diagnostics warrant sensitivity analyses and careful reporting.

## Interpretation and limitations

What this run establishes:
- A reproducible TCGA-trained, METABRIC-tested clinical-only Cox baseline ran successfully.
- The external METABRIC cohort was not used for model fitting or tuning.
- The model achieved Harrell C-index 0.616 (bootstrap 95% interval 0.592–0.641) on 1,454 complete-case external patients.
- Aggregate data-quality and exclusion counts are recorded.

What it does not establish:
- Clinical validation, clinical utility, prospective performance or treatment benefit.
- Molecular cross-platform validation.
- Prediction of immunotherapy response or any treatment response.
- Validation of tumour growth PDEs, PINNs, agent-based cells, GNN interactions, oxygen gradients or drug diffusion.
- A causal effect of age or stage.

Important threats:
- The broad-stage mapping is a pragmatic harmonization of TCGA AJCC pathologic stage and METABRIC TUMOR_STAGE; it needs source-level review.
- METABRIC complete-case selection excludes 525 expression patients for stage missingness.
- Stage IV is represented by only 10 external patients.
- The endpoint uses cBioPortal OS conventions; source-specific time origin and coding should be confirmed against the cohort documentation.
- This model uses clinical covariates only. Molecular data were not used.

## Next validation work

1. Review the exact TCGA and METABRIC stage-value dictionaries and document the broad-stage mapping, including ambiguous/unknown categories.
2. Add sensitivity analyses: age-only, stage-only, complete-case versus alternative missingness strategies, and stage grouping variants.
3. Assess calibration at prespecified time horizons and compare with a clinical baseline appropriate to each cohort.
4. Only after stage and endpoint harmonization is locked, design a molecular score with a fixed gene set and a prespecified cross-platform transformation.
5. Evaluate the molecular score in METABRIC without using METABRIC outcomes to select genes, thresholds or transformations.

## References

- cBioPortal FAQ: https://docs.cbioportal.org/user-guide/faq/
- cBioPortal clinical field conventions: https://github.com/cbioportal/cbioportal/blob/master/docs/File-Formats.md
- TCGA PanCancer Atlas: https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018
- METABRIC: https://www.cbioportal.org/study/summary?id=brca_metabric


## Stage and missingness sensitivity analysis (8 November 2026)

### Reproducibility

- Sensitivity workflow: https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/42455733551
- Workflow commit: `96407b83108b85ab81a5ceb1f99f0547bf5f322e`
- Conclusion: success
- Aggregate report artifact: `tcga-metabric-clinical-survival-baseline` (artifact ID `13607636020`)

The primary age + four-level stage model was left unchanged. Three additional analyses were run against the same external METABRIC endpoint.

### Results

| Sensitivity model | METABRIC patients | Harrell C-index | Bootstrap 95% interval |
|---|---:|---:|---:|
| Age only, stage-complete cases | 1,454 | 0.5914 | 0.5700–0.6174 |
| Age + stage I-IV (primary) | 1,454 | 0.6161 | 0.5923–0.6411 |
| Stage only, I reference | 1,454 | 0.4782 | 0.4565–0.4984 |
| Age + stage I-II vs III-IV | 1,454 | 0.6213 | 0.5966–0.6460 |
| Age only, no stage requirement | 1,979 | 0.5882 | 0.5683–0.6066 |

The no-stage-requirement age-only analysis retained 1,069 TCGA development patients (151 deaths) and 1,979 METABRIC external patients (1,143 deaths). This is a different, larger analysis population and should not be treated as a like-for-like comparison with the complete-case model.

### Interpretation

- The age-only C-index is similar when stage is not required (0.588) versus the stage-complete subset (0.591). This suggests that the large stage-related exclusion does not materially change age-only discrimination in this particular run, but it does not eliminate selection bias.
- Collapsing stage III and IV gives a slightly higher C-index (0.621) than the four-level primary model (0.616). Their bootstrap intervals overlap; this does not establish a reliable improvement.
- The stage-only model has an external C-index below 0.5 (0.478). This is a warning that stage coding, stage definitions, transportability or the prognostic relation may not align across cohorts. It must not be interpreted as evidence that advanced stage is protective. Raw stage category distributions and source definitions need to be audited before further interpretation.
- These sensitivity results remain exploratory clinical-only analyses. They do not validate a molecular signature, predict treatment benefit or validate the spatial digital twin.

### Next required check

A new aggregate-only cBioPortal audit has been added to report the raw value frequencies for TCGA `AJCC_PATHOLOGIC_TUMOR_STAGE` and METABRIC `TUMOR_STAGE`. This is intended to establish whether the two fields encode comparable categories before any new model is built. Do not infer harmonization from the shared field names or from the current broad-stage mapping alone.
