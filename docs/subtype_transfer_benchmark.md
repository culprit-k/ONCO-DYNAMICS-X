# Five-class intrinsic-subtype transfer benchmark

**Implementation:** `evaluate_intrinsic_subtype_transfer` in `src/onco_dynamics_x/data_engine/molecular_model.py`  
**Runner:** `scripts/run_pam50_molecular_validation.py`  
**Workflow:** `.github/workflows/pam50-molecular-validation.yml`  
**Status:** implementation and synthetic unit tests passed; public-cohort subtype metrics are pending the current GitHub Actions evaluation run.

## Question

Can a fixed 50-gene expression panel, transformed without using subtype labels, learn five existing breast-cancer subtype annotations in TCGA-BRCA and transfer those labels to the independent METABRIC cohort?

This is a supervised cross-cohort annotation-transfer benchmark. It is separate from the overall-survival model.

## Label definitions

- TCGA development target: patient-level `SUBTYPE` annotation, with values such as `BRCA_Basal`, `BRCA_Her2`, `BRCA_LumA`, `BRCA_LumB`, and `BRCA_Normal`.
- METABRIC external target: patient-level `CLAUDIN_SUBTYPE` annotation (PAM50 + claudin-low), with the five shared labels plus additional values such as `claudin-low` and `NC`.

Labels are canonicalised to the five shared classes:

```
Basal, Her2, LumA, LumB, Normal
```

METABRIC claudin-low and NC annotations are counted and explicitly excluded from this five-class evaluation rather than being forced into a PAM50 class. Patients without an annotation are not scored. The task uses the primary-tumour sample and survival-linked patient cohort selected by the main molecular pipeline.

## Features and leakage control

The classifier uses the 50-gene panel values after the same fixed, per-gene, cohort-wise mid-rank inverse-normal transform described in [the molecular methodology](pam50_molecular_validation.md).

- Expression transformation is completed without subtype labels or survival outcomes.
- The full unlabeled expression distribution of each cohort is used for rank transformation; this is transductive and not a frozen prospective transformation.
- TCGA labels are used only as supervised training targets.
- METABRIC subtype labels are held out and used only for external scoring.
- No METABRIC subtype labels are used for feature selection, scaling, model tuning or fitting.

## Classifier

The prespecified classifier is L2-regularised logistic regression with:
- (C=0.1);
- balanced class weights;
- LBFGS solver;
- maximum 3,000 iterations;
- no hyperparameter search.

The baseline is a constant classifier that predicts the most frequent shared subtype in the TCGA training labels.

## Metrics

External METABRIC outputs:
- accuracy;
- balanced accuracy;
- macro-F1;
- patient-bootstrap 95% intervals for those metrics;
- per-class precision, recall, F1 and support;
- five-by-five confusion matrix;
- raw annotation counts and explicit counts of excluded claudin-low/NC labels.

All values are aggregate-only. Patient IDs, sample IDs and row-level subtype labels are not written to the report.

## Interpretation limits

This benchmark uses existing cohort subtype labels as its target. Because the panel features overlap the historical intrinsic-subtype gene panel, high accuracy would primarily demonstrate that the fixed panel encodes the cohort annotations across these assays. It would **not** be independent biological discovery or clinical validation.

It does not reproduce official PAM50 centroids, commercial Prosigna, Risk-of-Recurrence scores or clinical thresholds. It does not predict treatment response, survival, recurrence or drug benefit. It does not establish that the PDE/agent-based digital twin correctly represents tumour biology.

## Acceptance gates after the external run

1. Check class support and annotation coverage in each cohort.
2. Compare accuracy, balanced accuracy and macro-F1 with the training-majority baseline.
3. Inspect per-class errors, especially confusion between LumA and LumB and between HER2-enriched and Basal.
4. Repeat with an assay-aware frozen transformation and probe-level METABRIC mapping before claiming prospective cross-platform performance.
5. Confirm on another independent cohort with compatible subtype annotation.
