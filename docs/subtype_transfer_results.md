# Five-class molecular subtype transfer — completed results

**Status:** completed successfully  
**Date:** 11 November 2026  
**Workflow:** [PAM50-panel molecular external validation run 42987068037](https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/42987068037)  
**Evaluated commit:** `05ce6675ff89b305acf27b75b2169713be9bc371`  
**Aggregate report:** `results/molecular_validation/pam50_tcga_metabric_summary.json`

## Task and cohort labels

A five-class classifier was trained on TCGA-BRCA patient-level `SUBTYPE` labels and evaluated once on METABRIC patient-level `CLAUDIN_SUBTYPE` labels. The five shared classes were Basal, Her2, LumA, LumB and Normal.

Training cohort:
- 968 patients with a survival-linked expression row and one of the five shared subtype labels.
- Basal 169, Her2 75, LumA 496, LumB 193, Normal 35.

External test cohort:
- 1,756 METABRIC patients with one of the five shared labels.
- Basal 209, Her2 224, LumA 700, LumB 475, Normal 148.
- Additional METABRIC annotations excluded from this five-class task: claudin-low 217, NC 6.

The task uses only the existing survival-linked, verified Primary Solid Tumor sample cohort from the main molecular workflow. Patients without an annotation are not scored. The five-class model does not attempt to classify claudin-low or NC.

## Model

- Features: all 50 PAM50 research-panel genes after within-cohort, outcome-blind rank-normalisation.
- Classifier: L2 logistic regression.
- Fixed (C=0.1), balanced class weights, LBFGS solver, maximum 3,000 iterations.
- No hyperparameter search and no METABRIC label use in training or tuning.
- Baseline: constant prediction of the most frequent subtype in the TCGA training labels (LumA).

## External METABRIC performance

| Metric | Result | Patient-bootstrap 95% CI |
|---|---:|---:|
| Accuracy | **81.0%** | 79.1–82.7% |
| Balanced accuracy | **74.7%** | 72.4–76.8% |
| Macro-F1 | **77.8%** | 75.3–79.9% |
| Majority-class baseline accuracy | 39.9% | — |

### Per-class performance

| Class | Precision | Recall | F1 | Test support |
|---|---:|---:|---:|---:|
| Basal | 98.2% | 78.9% | 87.5% | 209 |
| Her2 | 80.9% | 73.7% | 77.1% | 224 |
| LumA | 76.1% | 95.0% | 84.5% | 700 |
| LumB | 87.5% | 73.9% | 80.1% | 475 |
| Normal | 70.6% | 52.0% | 59.9% | 148 |

### Confusion matrix

Rows are the existing METABRIC label; columns are the model prediction.

| Actual  Predicted | Basal | Her2 | LumA | LumB | Normal |
|---|---:|---:|---:|---:|---:|
| Basal | 165 | 22 | 4 | 3 | 15 |
| Her2 | 1 | 165 | 23 | 32 | 3 |
| LumA | 0 | 6 | 665 | 15 | 14 |
| LumB | 0 | 9 | 115 | 351 | 0 |
| Normal | 2 | 2 | 67 | 0 | 77 |

The strongest performance was for Basal and LumA. Normal-like was the weakest class; 67 of 148 Normal-labelled cases were assigned LumA. LumB-to-LumA confusion was also common (115 cases). These patterns should be investigated in the context of assay differences and the study-specific annotation definitions.

## Interpretation

The five-class model transfers the existing study annotations substantially better than a constant-majority prediction on this external cohort. This is useful evidence that the fixed 50-gene panel retains subtype-associated information across the two gene-expression profiles.

However, the labels are existing study-provided annotations, and the panel features overlap the historical intrinsic-subtype gene inventory. This is an annotation-transfer benchmark, not independent biological discovery. The classifier can partly recover the algorithmic or biological structure already encoded in the labels; high performance does not validate an entirely new subtype taxonomy.

The result does not reproduce the official PAM50 centroids, commercial Prosigna assay, Risk-of-Recurrence score, or its clinical thresholds. It does not establish treatment response, recurrence prediction, clinical utility or individual-patient deployment.

## Limitations

- The feature transformation remains cohort-wise and transductive: it uses the full unlabelled expression distribution of each cohort.
- TCGA RNA-seq RSEM and METABRIC Illumina HT-12 v3 log2 intensities remain different platforms. The gene-level cBioPortal profiles were used; the original METABRIC probe-level annotation and duplicate-probe summarisation were not independently reconstructed.
- METABRIC claudin-low and NC categories are excluded from the five-class task.
- The task uses survival-linked patients with available subtype annotation; patients without OS follow-up or without subtype labels are not part of this classifier evaluation.
- Study annotations are not an independent histopathological ground truth.
- This task does not validate any spatial PDE, agent-based, immune-evasion, hypoxia or simulated-drug-response component.

## Next steps

1. Reconstruct the METABRIC probe-to-gene mapping and duplicate-probe summarisation from the original array data.
2. Replace the cohort-wise rank transform with a locked training-reference transform suitable for prospective single-sample application.
3. Evaluate the frozen classifier on another independent cohort with compatible annotation.
4. Analyse uncertainty and class-specific failure modes, especially Normal versus LumA and LumB versus LumA.
5. Keep the survival model and subtype classifier as separate tasks; do not imply subtype classification accuracy demonstrates prognostic accuracy or treatment response.
