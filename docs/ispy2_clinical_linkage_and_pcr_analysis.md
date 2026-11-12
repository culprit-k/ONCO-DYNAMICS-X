# I-SPY 2 feature–clinical linkage and exploratory pCR analysis

**Status: 12 November 2026.** The uploaded 985-person clinical workbook was matched against the 384-person multi-feature MRI workbook using the exact participant key. All 384 MRI rows matched one-to-one; no duplicates or unmatched MRI participants were found. The other 601 clinical rows have no corresponding feature row.

## Linkage

- Clinical workbook: 985 rows and 985 unique `Patient_ID` values.
- MRI feature workbook: 384 rows and 384 unique `CLINICAL-TRIAL-SUBJECT-ID` values.
- Exact one-to-one matches: **384/384 (100%)**.
- Duplicate IDs in either file: **0**.
- Matched participants missing Arm, HR, HER2, MP, pCR or age: **0**.

Join key: `CLINICAL-TRIAL-SUBJECT-ID = Patient_ID`. The 601 clinical-only records were not included in MRI analyses.

## Outcome and subgroup summaries

In the linked cohort, pCR code 1 appears in 113/384 (29.43%) and code 0 in 271/384 (70.57%). The original 384-participant MRI paper specifies binary pCR coding as 1 = pCR and 0 = non-pCR; this is the coding basis used here. Retain the variable-specific mapping when transferring to another release or cohort.

| HR | HER2 | n | pCR code 1 | Rate |
|---:|---:|---:|---:|---:|
| 0 | 0 | 132 | 50 | 37.9% |
| 0 | 1 | 30 | 20 | 66.7% |
| 1 | 0 | 162 | 24 | 14.8% |
| 1 | 1 | 60 | 19 | 31.7% |

The cohort contains 13 treatment-arm labels. Individual arm sizes range from 9 to 65, so arm-specific pCR rates are descriptive and should not be interpreted as comparative efficacy or causal treatment effects.

## FTV-candidate change by pCR code

| Visit | Median change, code 0 | Median change, code 1 |
|---|---:|---:|
| T1 | -39.14% | -62.31% |
| T2 | -82.91% | -93.25% |
| T3 | -93.53% | -96.56% |

Mann–Whitney comparisons of log volume ratios yielded p=1.99e-9 (T1), p=8.13e-11 (T2), and p=8.79e-6 (T3). These are exploratory unadjusted associations; no multiplicity correction was applied.

## Logistic association

Unadjusted odds ratios are per one-SD increase in `log(volume_visit / volume_T0)`. A larger log ratio means less volume reduction.

| Visit | OR | 95% CI | p |
|---|---:|---:|---:|
| T1 | 0.532 | 0.417–0.678 | 3.44e-7 |
| T2 | 0.471 | 0.365–0.606 | 5.34e-9 |
| T3 | 0.596 | 0.471–0.754 | 1.63e-5 |

A larger volume ratio (less reduction) was associated with lower odds of pCR code 1. This is an association, not a causal treatment effect. A secondary T1 model adjusted for age, HR and HER2 gave OR 0.525 (95% CI 0.405–0.680) for the standardized log ratio; the T2 adjusted OR was 0.429 (0.326–0.565). Treatment arm was not included in those coefficient models.

## Internal cross-validation

Repeated stratified 5-fold CV, 10 repeats (50 folds), gave:

| Model | Mean ROC AUC | Fold SD |
|---|---:|---:|
| T1 log ratio + age + HR + HER2 | 0.753 | 0.050 |
| Same + treatment arm | 0.771 | 0.055 |
| T2 log ratio + age + HR + HER2 | 0.771 | 0.039 |
| Same + treatment arm | 0.790 | 0.043 |

These are internal estimates in the same 384-person cohort, not external validation or clinical utility evidence. Fold SD is descriptive, not a confidence interval. T3 is post-neoadjuvant/pre-surgery and is not used as an early prediction feature.

## Interpretation and remaining gates

This analysis now supports cohort-level association work between serial MRI features and pCR. It does **not** provide elapsed-day kinetics because actual acquisition dates are absent. It also does not establish the volume field's exact physical units or equivalence to trial FTV, nor does it make FTV a direct measurement of viable cancer-cell number.

Before mechanistic calibration or publication:
1. Confirm pCR coding, HR/HER2/MP definitions and volume units against the official data dictionary.
2. Obtain acquisition dates and drug exposure timelines if available.
3. Pre-specify a parsimonious prediction model and treatment-arm adjustment strategy; use locked validation and calibration assessment.
4. Only after verified units and time intervals, estimate treatment-exposed phenomenological kinetics. Do not label them cell birth/death rates without a validated observation model.
5. Keep the participant-level linked workbook private. This repository contains aggregate results and code only.

## Follow-on response model

A first patient-separated model comparison has been completed using the same linked cohort. The 80/20 stratified split was made on unique participant identifiers (307 train, 77 test; zero patient overlap). A linear response reference and a quadratic nonlinear response surface were evaluated using T1/T2 log-volume ratios, age, HR, HER2 and treatment arm. The pCR classifier did not use T3; a separate model forecast T3 response without using pCR. Test pCR AUC was 0.726 for the linear reference and 0.739 for the quadratic model; the paired difference was +0.013 (95% bootstrap CI -0.018 to +0.050). T3 forecast RMSE was 0.997 versus 1.001. These results do not establish a reliable nonlinear-model advantage.

The coupled reaction-diffusion model was not numerically fitted. Whole-lesion volume and scalar morphology summaries do not identify spatial diffusion under no-flux boundaries, and actual acquisition times, treatment exposure timelines, registered spatial fields and a validated MRI observation operator are not available in the linked workbooks. See [the full response-model comparison and mathematical identifiability argument](ispy2_pcr_response_model_comparison.md). No elapsed-day rates or biological PDE parameters are reported.

## Data provenance

- TCIA identifies the Multi-feature MRI NACT workbook as the 384-participant dataset associated with the I-SPY 2 MRI response paper and separately lists the 985-participant clinical release: [TCIA I-SPY 2 collection](https://www.cancerimagingarchive.net/collection/ispy2/).
- The original multi-feature MRI study describes the imaging response analysis and its clinical endpoint context: [Li et al., npj Breast Cancer (2020)](https://doi.org/10.1038/s41523-020-00203-7).
- TCIA's collection page lists CC BY 4.0 for the clinical and multi-feature workbook resources. Preserve the collection citation and attribution in any redistributed analysis.

## Reproduction

Install the optional science and survival dependencies, then run:

```bash
pip install -e ".[science,survival]"
python scripts/analyze_ispy2_clinical_linkage.py \\
  --features /private/Multi-feature-MRI-NACT-Data.xlsx \\
  --clinical /private/ISPY2-Imaging-Cohort-1-Clinical-Data.xlsx \\
  --output /private/ispy2_clinical_linkage_and_pcr_aggregate.json
```

The script refuses duplicate or unmatched feature IDs and writes aggregate JSON only. It never writes participant-level linked rows.
