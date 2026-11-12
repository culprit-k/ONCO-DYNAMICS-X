# I-SPY 2 MRI workbook: initial aggregate analysis

**Analysis date:** 12 November 2026  
**Dataset:** Multi-feature MRI NACT Data.xlsx  
**Privacy:** Aggregate-only output; participant identifiers and row-level records are excluded.

## Dataset audit

- 384 rows and 29 columns.
- 384 unique subject identifiers; no duplicate subject rows.
- No missing cells in the workbook.
- All 384 subjects have four positive volume measurements.
- The three supplied FTV percentage-change fields agree with the corresponding volume-derived percentages to floating-point precision.

## Longitudinal summary

| Visit | Geometric mean | Bootstrap 95% CI |
|---|---:|---:|
| T0 | 15.744 | 14.181–17.535 |
| T1 | 7.507 | 6.597–8.525 |
| T2 | 1.929 | 1.619–2.275 |
| T3 | 0.836 | 0.712–0.987 |

The geometric-mean T3/T0 ratio was **0.0531** (paired bootstrap 95% CI 0.0462–0.0608), corresponding to a geometric-mean reduction of **94.69%**. The median individual percentage change was -94.58%; 378/384 subjects had a lower value at T3 than T0, while 6/384 had a higher value.

## Bayesian endpoint estimate

A conjugate reference-prior Normal model on the 384 subject-level log(T3/T0) ratios gives a posterior mean log ratio of **-2.935** and a 95% credible interval for the geometric ratio of **0.0462–0.0611**. This is an endpoint contrast, not a time-normalized biological rate.

## What cannot yet be estimated

The workbook has no actual acquisition dates, treatment arm/regimen, HR/HER2 subtype or pCR outcome. The four MRI visits also span different treatment stages and unequal intervals. Therefore this analysis does **not** estimate treatment-adjusted effects, per-day kinetic rates, drug response parameters, PDE growth/death coefficients, or predictive performance.

The workbook's column headers do not state the physical volume unit. Confirm the exact semantic mapping and units of VOLUME_TUM_BLU_V10–V40 against the source data dictionary before assigning physical units or integrating the values into a tumour-state variable.

## Required next inputs

1. The matching clinical workbook/data extract with regimen, subtype, pCR and documented linkage key.
2. Actual acquisition dates or elapsed days for each MRI visit.
3. Data dictionary confirmation of volume field definition and units.

The trial protocol describes T0 as pretreatment, T1 as early treatment after three cycles, T2 after twelve cycles/between regimens, and T3 post-AC/pre-surgery. These intervals are not equal. The clinical paper describes a 384-person complete-case analysis and reports the same four MRI features; it does not make this workbook a randomized treatment-arm dataset for causal comparisons.

All participant-level linkage and source data must remain in approved local storage; only aggregate reports belong in the public repository.
