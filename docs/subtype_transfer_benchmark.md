# Five-class intrinsic-subtype transfer benchmark

**Implementation:** `evaluate_intrinsic_subtype_transfer` in `src/onco_dynamics_x/data_engine/molecular_model.py`  
**Runner:** `scripts/run_pam50_molecular_validation.py`  
**Workflow:** `.github/workflows/pam50-molecular-validation.yml`  
**Status:** completed successfully on public TCGA-BRCA and METABRIC cohorts; full tests and direct Ruff checks passed.

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

## Completed external results

**Workflow:** [run 42987068037](https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/42987068037)

- Training: 968 TCGA patients with one of the five shared labels.
- External test: 1,756 METABRIC patients with one of the five shared labels.
- METABRIC extra labels excluded: claudin-low 217; NC 6.
- Accuracy: **81.0%** (bootstrap 95% CI 79.1–82.7%).
- Balanced accuracy: **74.7%** (95% CI 72.4–76.8%).
- Macro-F1: **77.8%** (95% CI 75.3–79.9%).
- Majority-class baseline accuracy: 39.9%.

The weakest class was Normal-like (F1 0.60), driven largely by confusion with LumA. LumB-to-LumA confusion was also common. See [full result and confusion matrix](subtype_transfer_results.md) and the [aggregate JSON report](../results/molecular_validation/pam50_tcga_metabric_summary.json).

## Interpretation limits

This benchmark uses existing cohort subtype labels as its target. Because the panel features overlap the historical intrinsic-subtype gene panel, high accuracy would primarily demonstrate that the fixed panel encodes the cohort annotations across these assays. It would **not** be independent biological discovery or clinical validation.

It does not reproduce official PAM50 centroids, commercial Prosigna, Risk-of-Recurrence scores or clinical thresholds. It does not predict treatment response, survival, recurrence or drug benefit. It does not establish that the PDE/agent-based digital twin correctly represents tumour biology.

## Next validation gates

1. Verify the original TCGA and METABRIC subtype annotation definitions and versions, including how Normal-like, claudin-low and NC were assigned.
2. Rebuild the METABRIC probe-to-gene mapping and duplicate-probe summarisation from the original Illumina HT-12 v3 annotation rather than relying only on the cBioPortal gene-level profile.
3. Replace the cohort-wise transductive rank transformation with a frozen reference transform learned in training data and usable on one new sample at a time.
4. Confirm transfer on another independent cohort with a compatible assay and comparable subtype annotations.
5. If the goal becomes clinical use, define an intended use and conduct a separate prospective clinical-validation and utility study.
