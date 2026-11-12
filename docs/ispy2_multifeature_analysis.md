# I-SPY 2 multi-feature longitudinal analysis

**Status: 12 November 2026.** The local 384-record Multi-feature MRI NACT Data workbook was analysed with the aggregate-only runner `scripts/analyze_ispy2_multifeature.py`. The raw workbook, participant identifiers and row-level measurements are not part of this repository.

## Reproduce locally

```bash
pip install -e ".[science,dev]"
python scripts/analyze_ispy2_multifeature.py \
  --input /path/to/Multi-feature-MRI-NACT-Data.xlsx \
  --output results/ispy2/ispy2_multifeature_aggregate.json \
  --bootstrap 2000 \
  --seed 20261112
pytest tests/test_ispy2_multifeature_analysis.py
```

The script validates uniqueness of participant keys, rejects missing/non-finite or negative measurements, verifies every supplied percentage-change field against the corresponding baseline and follow-up measurements, and exports aggregate summaries only. It uses Friedman repeated-measures tests for the four ordered visits, paired Wilcoxon tests for each follow-up versus baseline, Holm adjustment across the four omnibus tests and separately across the twelve planned paired contrasts, percentile bootstrap intervals for median changes, and exploratory Spearman associations between endpoint changes.

## Key descriptive results

All 384 participants had four complete observations for each of the four feature families. All twelve supplied percentage-change fields reproduced from the raw measurement columns to a maximum absolute error below 1e-6 percentage points.

| Feature | T3 vs T0 median change | Participants lower at T3 | Friedman p-value |
|---|---:|---:|---:|
| Workbook volume field (FTV candidate) | -94.58% | 378/384 | 4.23 × 10⁻¹⁹¹ |
| Longest diameter | -67.05% | 365/384 | 5.80 × 10⁻¹⁶⁵ |
| Sphericity | +25.55% | 274/384 | 4.32 × 10⁻²⁹ |
| Background parenchymal enhancement | -28.95% | 294/384 | 2.86 × 10⁻⁴⁸ |

The Friedman p-values test for any within-subject difference across the four visit stages; they do not establish a monotonic trend, causal effect or clinically meaningful change. For the T3-vs-T0 paired tests, p-values were adjusted jointly across the twelve feature-by-visit contrasts using Holm's method. The aggregate report contains the corresponding confidence intervals, counts and adjusted values.

The volume-field percentage change was moderately associated with sphericity change (Spearman ρ = -0.509; bootstrap 95% interval -0.581 to -0.430), longest-diameter change (ρ = 0.365; 0.272 to 0.457), and more weakly with background parenchymal enhancement change (ρ = 0.245; 0.151 to 0.337). These are exploratory cross-feature associations, not mechanistic couplings or independent predictors.

## Measurement and inference boundaries

- TCIA describes the 384-person workbook as containing functional tumour volume, longest diameter, sphericity and contralateral background parenchymal enhancement. The volume columns in the workbook are named `VOLUME_TUM_BLU_V10` through `V40`; the script deliberately calls them the **workbook volume field / FTV candidate** until their exact mapping and unit are confirmed from the workbook-specific documentation.
- The separate TCIA DICOM data dictionary defines DICOM FTV volume in cc. This does **not**, by itself, prove that the workbook columns have identical semantics, units, segmentation platform or implementation. TCIA also warns that the UCSF-derived DICOM FTV can differ from the I-SPY 2 trial's Hologic/AEGIS FTV. Do not silently equate these measurements.
- T0 is pre-treatment; T1 is early treatment after about three weeks; T2 is mid-treatment between regimens; T3 is post-treatment before surgery. They are treatment-stage labels, not equal-duration time intervals. The present workbook does not contain actual scan dates, drug-exposure timing, treatment arm, subtype or pCR.
- The cohort is treatment-exposed. The observed decline is not an untreated growth rate, randomized treatment-effect estimate, drug-killing parameter, EC50, or direct viable-cancer-cell measurement.
- The current analysis does not estimate any ONCO-DYNAMICS-X PDE coefficient. In particular, whole-lesion enhancement-volume change cannot identify spatial diffusion, immune killing, hypoxia, vascular supply or fibroblast-transition parameters.

## Official sources

- [TCIA I-SPY 2 collection, cohort and data access](https://www.cancerimagingarchive.net/collection/ispy2/)
- [Li et al. (2020), multi-feature MRI response analysis](https://doi.org/10.1038/s41523-020-00203-7)
- [TCIA ACRIN 6698 / I-SPY 2 DWI and DCE MRI data dictionary](https://www.cancerimagingarchive.net/wp-content/uploads/ACRIN-6698-ISPY2-DWI-and-DCE-MRI-Data-Descriptions_20210520.pdf)

The current TCIA I-SPY 2 record lists CC BY 4.0 for the 384-person workbook and the separate 985-case clinical dataset. Cite the dataset and publication when using them.
