# Phase 3 molecular external validation — results

**Analysis:** fixed PAM50 50-gene panel, TCGA-BRCA development / METABRIC external evaluation  
**Status:** completed exploratory external evaluation  
**Date:** 11 November 2026  
**Successful verified-sample GitHub Actions run:** [42985971698](https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/42985971698)  
**Evaluated commit:** `9d14003f04602fe40e10a2c1070c2b7dffaa7824)  
**Aggregate result:** `results/molecular_validation/pam50_tcga_metabric_summary.json`

## Question

Does a predeclared breast-cancer expression panel contribute measurable overall-survival discrimination beyond age at diagnosis when a model fitted in TCGA-BRCA is transported to the independent METABRIC cohort?

This is a retrospective, exploratory prognostic-association analysis. It is not a clinical tool, a treatment-response predictor, or a validation of the ONCO-DYNAMICS-X spatial/mechanistic digital twin.

## Cohorts and data quality

| Cohort | Expression profile samples | Patients in survival analysis | Death events | Censored |
|---|---:|---:|---:|---:|
| TCGA-BRCA PanCancer Atlas | 1,082 | 1,069 | 151 | 918 |
| METABRIC | 1,980 | 1,979 | 1,143 | 836 |

All 50 genes in the fixed research panel had finite values in both expression profiles:
- TCGA: 54,100 / 54,100 expected sample–gene measurements.
- METABRIC: 99,000 / 99,000 expected sample–gene measurements.
- No duplicate sample–gene records, non-finite values, or gene-value missingness were observed.

Historical panel identifiers were resolved against current cBioPortal gene records: CDCA1 to NUF2, KNTC2 to NDC80, and ORC6L to ORC6.

### Sample selection verified

The final successful run used the cBioPortal DETAILED sample projection and parsed the official Sample fields (sample ID, sample type and nested patient ID). Every sample in both selected expression profiles was returned as `Primary Solid Tumor`:

- TCGA-BRCA: 1,082 / 1,082 primary solid tumour samples.
- METABRIC: 1,980 / 1,980 primary solid tumour samples.
- Samples without patient linkage: zero in both cohorts.
- Duplicate primary samples excluded: zero in both cohorts.

The external model results below were reproduced on this verified primary-solid-tumour sample set. See [successful corrected workflow run 42985971698](https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/42985971698) and the [aggregate JSON report](../results/molecular_validation/pam50_tcga_metabric_summary.json).

Subtype annotations were available in the two cohorts (TCGA `SUBTYPE`, METABRIC `CLAUDIN_SUBTYPE`). They were audited only as aggregate counts; no subtype labels were used for expression transformation, PCA, model fitting, tuning, or evaluation.

## Molecular feature construction

The fixed input inventory was the published 50-gene PAM50 research panel. The project did not use official PAM50 centroids, a commercial Prosigna implementation, Risk-of-Recurrence scores, or clinical cut-points.

1. Gene symbols were resolved through the cBioPortal gene catalogue.
2. Each gene was rank-normalised independently within each cohort using mid-ranks and an inverse-normal transform. This is outcome-blind but uses the full unlabelled expression distribution within each cohort.
3. Five principal components were fitted using TCGA only.
4. PCA loadings and component score scaling were frozen from TCGA and applied unchanged to METABRIC.
5. No survival outcomes or subtype annotations from METABRIC were used to select genes, construct features, or tune the model.

The five TCGA-fitted components explained 68.19% of variance in the rank-normalised TCGA panel matrix. PC1 had the strongest positive loadings on KIF2C, CDC20, CEP55, NDC80, MELK and UBE2C, and negative loadings on MLPH, NAT1, PGR, ESR1, FOXA1 and MAPT. This is a descriptive statistical axis consistent with a proliferation-versus-luminal-expression contrast in this dataset; it is not a causal mechanism or an intrinsic-subtype assignment.

The rank transformation is **transductive**: the full unlabelled METABRIC expression distribution is used to calculate METABRIC gene ranks. It is not a single-patient-at-a-time transform and is not ready for prospective individual-sample deployment.

## Models

Both models were fitted in TCGA on the same patient set:
- Age-only penalised Cox proportional-hazards reference.
- Age plus five PCA scores from the fixed expression panel, using an L2 Cox penalizer of 0.5.

Age standardisation and PCA component scaling were estimated in TCGA and applied unchanged to METABRIC. Stage was excluded because the TCGA patient-level AJCC pathologic stage and the METABRIC sample-level generic tumour-stage field have not been established as semantically equivalent.

## External performance

### Harrell concordance index

| METABRIC model | Harrell C-index | Patient-bootstrap 95% CI |
|---|---:|---:|
| Age only | 0.5882 | 0.5683–0.6065 |
| Age + five molecular PCs | **0.5992** | 0.5798–0.6201 |

The absolute difference was **+0.0111** in favour of the molecular model. The paired patient-bootstrap 95% confidence interval for the difference was **−0.0019 to +0.0241**, which includes zero. The analysis therefore does **not** demonstrate a statistically conclusive improvement over age alone.

### Time-dependent external performance

The report uses cumulative/dynamic time-dependent AUC and inverse-probability-of-censoring weighted Brier scores. The censoring distribution is estimated marginally within METABRIC for evaluation correction only.

| Horizon | Model | IPCW AUC (95% CI) | IPCW Brier |
|---|---|---:|---:|
| 12 months | Age only | 0.5783 (0.4602–0.6844) | 0.01505 |
| 12 months | Age + molecular PCs | 0.6233 (0.4995–0.7335) | 0.01503 |
| 36 months | Age only | 0.5151 (0.4740–0.5613) | 0.10128 |
| 36 months | Age + molecular PCs | 0.5440 (0.4967–0.5865) | 0.10099 |
| 60 months | Age only | 0.5395 (0.5037–0.5719) | 0.17203 |
| 60 months | Age + molecular PCs | 0.5798 (0.5469–0.6100) | 0.17085 |

At 12 and 36 months, the molecular model's AUC intervals include 0.5. At 60 months, its interval is above 0.5, but this is a single retrospective horizon and does not establish clinical utility or generalisability. The Brier-score differences are small.

At 60 months, the molecular model's mean predicted survival was 0.8104 versus an IPCW observed estimate of 0.7801; this corresponds to average survival overprediction of about 3.0 percentage points.

### Model diagnostics

The training-cohort rank-based proportional-hazards tests were non-significant for age and the five PCA scores. A non-significant diagnostic does not prove that proportional hazards holds.

No individual PCA coefficient was nominally significant at the 0.05 level in the fitted TCGA Cox model. The five PCs should therefore not be presented as independent, confirmed biological risk factors.

## Interpretation

The fixed panel is measurable in both public profiles, and the leakage-conscious feature pipeline executed. The age-plus-PC model showed a small apparent increase in external discrimination, but the paired uncertainty interval includes no improvement. Time-dependent results are mixed, and five-year survival calibration remains imperfect.

The appropriate conclusion is:

**The fixed expression panel is feasible to investigate, but its added prognostic value over age alone remains unconfirmed.**

This result should be treated as an exploratory baseline for further assay-level and cohort-definition work, not as a validated molecular signature.

## Limitations and next acceptance gates

1. **Sample type:** resolve specimen type in both cohorts and repeat the evaluation using verified primary-tumour samples only.
2. **Platform mapping:** inspect the original METABRIC Illumina HT-12 v3 probe annotation, multiple-probe mapping, duplicate summarisation, and annotation version. The present analysis uses cBioPortal's gene-level profile rather than rebuilding probe-level expression.
3. **Harmonisation:** compare a predeclared, reproducible assay-aware transformation with the current cohort-wise rank transform. Retain the current transductive limitation until a prospective transform exists.
4. **External validation:** repeat on a further independent cohort with a compatible, documented assay and endpoint.
5. **Clinical endpoint:** verify endpoint time origin and definitions against the original cohort documentation; assess informative or covariate-dependent censoring.
6. **Model utility:** compare calibration, time-dependent discrimination, prediction error, and clinical net benefit only after a defined intended-use question and independent validation design.
7. **Biological mechanism:** bulk expression and survival data do not validate cell positions, immune-cell contacts, oxygen gradients, spatial immune escape, drug penetration, or PDE/agent-based trajectories. Those require matched spatial/pathology data.

## Data governance

The public repository contains only aggregate counts, metrics, model summaries and PCA loading summaries. No patient IDs, sample IDs, raw expression values, or row-level clinical records were written or uploaded. Controlled/restricted data must remain in approved storage.

## References

- Parker JS et al. *Supervised Risk Predictor of Breast Cancer Based on Intrinsic Subtypes*. Journal of Clinical Oncology. 2009;27:1160–1167. https://doi.org/10.1200/JCO.2008.18.1370
- `genefu` PAM50 gene mapping (gene inventory source only; official centroids are not copied): https://github.com/bhklab/genefu/blob/master/inst/extdata/pam50_model.csv
- cBioPortal API and clients: https://github.com/cBioPortal/cbioportal/blob/master/docs/web-API-and-Clients.md
- cBioPortal clinical and molecular data conventions: https://docs.cbioportal.org/file-formats/
