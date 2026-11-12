# I-SPY 2 workbook audit and longitudinal calibration handoff

**Status (12 November 2026):** the 384-participant Multi-feature MRI NACT workbook is locally audited and analysed. A privacy-preserving structural audit, aggregate volume analysis and four-feature repeated-measures analysis are complete. The separate 985-case clinical workbook has not yet been locally obtained or linked. No treatment-adjusted effect, elapsed-day kinetic estimate or patient-derived PDE coefficient is claimed.

## Public data source

The TCIA I-SPY 2 collection currently lists:
- a 384-participant **Multi-feature MRI NACT Data** workbook (about 130 kB) containing DCE-derived FTV, longest diameter, sphericity and background parenchymal enhancement;
- a separate clinical spreadsheet for the 985-case I-SPY 2 Imaging Cohort 1 (719 I-SPY 2 plus 266 ACRIN-6698 cases);
- the 719-case DCE collection and a combined 985-case imaging manifest.

The 384-person feature workbook and 985-case clinical workbook are different resources. Do not assume they have identical identifiers, participants, visit coverage or treatment fields. Establish an explicit, approved linkage using the data dictionary and documented join keys before combining them. Do not upload participant-level workbooks or linkage tables to this public repository.

Source: [TCIA I-SPY 2 collection and data access](https://www.cancerimagingarchive.net/collection/ispy2/). The collection lists CC BY 4.0; retain the collection DOI/citation and recheck current access terms when downloading.

## Completed local analysis

The feature workbook contains 384 unique subject rows and 29 columns. All participants have four complete observations for each of the four feature families. The twelve supplied percentage-change fields reproduce their corresponding baseline/follow-up formulas with maximum absolute error below 1e-6 percentage points.

Median T3-versus-T0 changes were:
- workbook volume field (FTV candidate): -94.58%; 378 lower and 6 higher;
- longest diameter: -67.05%; 365 lower, 15 higher and 4 unchanged;
- sphericity: +25.55%; 274 higher and 110 lower;
- background parenchymal enhancement: -28.95%; 294 lower and 90 higher.

Friedman tests found within-subject differences across the four ordered visit stages for all four features. These are descriptive stage contrasts. They do not imply equal elapsed intervals, a causal treatment effect, or a biological PDE rate. See [the detailed multi-feature analysis](ispy2_multifeature_analysis.md) and [aggregate JSON](../results/ispy2/ispy2_multifeature_aggregate.json).

## Measurement and visit limitations

1. **Volume definition and units:** the workbook columns are named VOLUME_TUM_BLU_V10 through V40. The headers alone do not establish physical units or equivalence to DICOM-derived FTV. The separate TCIA DICOM dictionary defines DICOM FTV volume in cc, but this does not prove the workbook columns use the same implementation. TCIA also notes that UCSF-derived DICOM FTV can differ from the trial's Hologic/AEGIS measurement. Confirm the workbook-specific mapping before assigning units or using the existing cc-based fitter.
2. **Visit timing:** T0 is pretreatment; T1 is early treatment after about three weeks of regimen 1; T2 is between paclitaxel and AC; T3 is after AC and before surgery. These are not equally spaced days. Actual MRI acquisition dates are absent from the feature workbook.
3. **Treatment context:** the workbook lacks treatment arm/regimen, drug exposure timing, subtype, pCR and actual scan dates. The cohort is treatment-exposed, so the aggregate decrease cannot be treated as untreated tumour growth or a causal treatment effect.
4. **Measurement interpretation:** functional enhancement volume is not a direct measurement of viable cancer-cell count. Changes may include response-related enhancement changes and measurement/segmentation effects.

## Run the structural audit locally

Install the science dependencies and run the audit against each downloaded workbook separately:

```bash
pip install -e ".[science]"
python scripts/audit_ispy2_workbook.py \
  --input /path/to/Multi-feature-MRI-NACT-Data.xlsx \
  --output /path/to/private/ispy2_feature_workbook_audit.json

python scripts/audit_ispy2_workbook.py \
  --input /path/to/I-SPY2-clinical.xlsx \
  --output /path/to/private/ispy2_clinical_workbook_audit.json
```

The report includes sheet dimensions, column names, candidate semantic roles, data types, missingness and unique-value counts. It does not export cell values, row-level observations or participant identifiers. Column-role matching is heuristic and must be manually reviewed against the workbook's documentation.

## Required audit decisions before fitting

1. **Workbook identity:** record exact source title, version/download date, license and data citation locally.
2. **Visit semantics:** map every scan/visit field to actual acquisition day or date. Nominal labels such as T0/T1/T2/T3 are not elapsed time.
3. **FTV units and definition:** verify the workbook's own units and method. The DICOM-derived UCSF FTV and the trial's primary commercial-platform FTV are not silently interchangeable.
4. **Participant linkage:** verify the documented key, one-to-one or one-to-many cardinality, duplicate rows, unmatched cases and overlap between the 384-feature and 985-clinical workbooks. Keep the mapping private.
5. **Treatment exposure:** align regimen/arm, start and stop dates, drug sequence and MRI acquisition dates. A regimen label alone does not establish exposure at a particular scan.
6. **Subtype and response:** verify ER/PR/HER2 definitions, assay dates, pCR definition and whether the measurement is baseline or post-treatment.
7. **Temporal eligibility:** quantify participants with at least three positive FTV observations and actual acquisition days. Zero/below-detection FTV is censored, not zero-volume evidence.
8. **Holdout:** split by participant, not row. Predeclare an early-visit calibration window and a later-visit temporal holdout before model fitting.

## Estimand and modelling boundary

The initial estimand is the **treatment-exposed net log-FTV change per day**:

```
log(FTV_i(t)) = alpha_i + k_net,i * elapsed_days + error_i(t)
```

The per-series slope is descriptive of the observed treatment context. It does not separate intrinsic proliferation from therapy-induced death, and it is not a drug EC50, drug-killing coefficient, or the PDE `cancer_growth` parameter.

The existing runner `scripts/fit_treatment_net_volume.py` accepts a local long-format CSV after the workbook audit and manual mapping. Its default volume field is `ftv_cc`; only use that mapping after verifying the workbook's actual unit. Keep raw workbooks, participant linkage and patient-level outputs in approved local/institutional storage.

## Acceptance gate for the next modelling stage

Do not start cohort-level fitting until the local audit confirms:
- interpretable visit/date columns and verified units;
- an auditable participant linkage and duplicate policy;
- aligned treatment exposure metadata;
- enough eligible serial positive FTV measurements;
- a participant-level temporal holdout that does not leak later visits into calibration.

If these conditions fail, report the exact data limitation and revise the estimand rather than imputing unsupported dates, inventing visit intervals or fitting mechanistic PDE coefficients.
