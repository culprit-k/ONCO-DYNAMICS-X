# TCGA–METABRIC stage harmonization audit

## Status

<badge>Audit completed; semantic harmonization unresolved</badge>

The cBioPortal audit workflows completed successfully and emitted aggregate-only stage category counts. This is a metadata audit, not proof that the clinical fields are comparable.

## Raw fields and observed categories

### TCGA-BRCA PanCancer Atlas

- Field: `AJCC_PATHOLOGIC_TUMOR_STAGE`
- Level: patient
- Display name: Neoplasm Disease Stage American Joint Committee on Cancer Code
- Values:
  - Stage I: 89
  - Stage IA: 86
  - Stage IB: 6
  - Stage II: 6
  - Stage IIA: 355
  - Stage IIB: 255
  - Stage III: 2
  - Stage IIIA: 155
  - Stage IIIB: 28
  - Stage IIIC: 64
  - Stage IV: 19
  - Stage X: 14

### METABRIC

- Field: `TUMOR_STAGE`
- Level: sample
- Display name: Tumor Stage
- Values:
  - 0: 24
  - 1: 630
  - 2: 979
  - 3: 144
  - 4: 11

Counts above are all populated records for the field in each study audit, not only the final survival-analysis complete cases.

## Why this matters

TCGA's field explicitly identifies AJCC pathologic stage. METABRIC's field is a generic sample-level Tumor Stage attribute. Although the numeric values appear ordinal, the audit does not establish:
- which staging system METABRIC used;
- whether it is clinical or pathological stage;
- the timepoint at which it was recorded;
- what category 0 means;
- whether its categories are semantically equivalent to TCGA AJCC pathologic stages.

cBioPortal stores study-defined patient and sample clinical attributes. The portal's clinical data model does not make custom fields comparable merely because their labels look similar.

## Link to model sensitivity result

The age + stage model reached METABRIC Harrell C-index 0.616; collapsing stages I-II versus III-IV reached 0.621; the stage-only model reached 0.478 (bootstrap 95% interval 0.457-0.498). The below-0.5 stage-only result is a warning. It must not be interpreted as evidence that higher stage is protective.

## Required next action

Find and review the original METABRIC clinical data dictionary and publication's stage definition. Confirm category 0 and the stage assessment timepoint. If equivalence to TCGA AJCC pathologic stage cannot be demonstrated, stop stage-based cross-cohort transport and retain age-only as the baseline comparator.

## Audit provenance

- METABRIC: https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/42455802465
- TCGA: https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/42455802428
- Repository commit: `c5aa4883052d48b845a426084a533e7c8fbae932`
