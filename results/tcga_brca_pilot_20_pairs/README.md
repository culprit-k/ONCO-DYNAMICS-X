# TCGA-BRCA pilot result

This directory stores the first successful real-data expression baseline from ONCO-DYNAMICS X.

- Source workflow run: [42452095954](https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/42452095954)
- Source code commit: `0286ee3ba7f0ebad04cbb031a79d7e80470d0a6a`
- Cohort: TCGA-BRCA
- Subset: 20 paired patients, 40 open-access STAR-count expression files
- Matrix: 40 samples by the common GDC gene set, using `tpm_unstranded`
- Evaluation: patient-grouped 75/25 holdout and five-fold StratifiedGroupKFold
- Result: holdout balanced accuracy 1.00 and ROC-AUC 1.00; grouped CV balanced accuracy and ROC-AUC 1.00

The result is a small technical specimen-label benchmark, not clinical prediction. It must not be interpreted as diagnosis, prognosis, treatment response, causal discovery or validation of the spatial digital twin.

Raw expression files and the expression matrix are deliberately not stored in this repository. The run artifact retained the selected file IDs, sample metadata and checksums for 14 days. A larger 48-patient paired evaluation has also been launched; its results are separate and must replace this pilot only after successful completion and review.
