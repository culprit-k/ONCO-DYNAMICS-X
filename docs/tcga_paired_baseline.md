# TCGA-BRCA paired expression baseline

## Research question for this first real-data benchmark

Can a simple, leakage-controlled expression model distinguish the known GDC specimen labels **Primary Tumor** and **Solid Tissue Normal** in a small set of patients who contributed both specimen types?

This is a technical data-pipeline and modelling benchmark. It is not a clinical diagnostic model, does not predict future cancer in healthy people, does not predict survival or treatment response, and does not validate the spatial PDE/ABM digital twin.

## Cohort selection

The GitHub Actions workflow queries the NCI Genomic Data Commons for:
- Project: `TCGA-BRCA`
- Data type: `Gene Expression Quantification`
- Workflow: `STAR - Counts`
- Access: `open`
- Sample classes: `Primary Tumor` and `Solid Tissue Normal`

It queries up to 500 metadata records per class, then selects one expression file per patient for each specimen class among patients present in both sets. The selected subset selects up to 48 paired patients (up to 96 samples). The initial metadata audit found 499 primary-tumour sample records and 113 normal sample records, representing 48 patients shared across the two queried sets. The workflow therefore selects up to 48 matched patient pairs from the current query window. The exact selected patient/sample IDs, GDC file UUIDs, file names, sizes and checksums are recorded in the workflow artifact. No expression matrices or raw expression files are committed to Git.

A patient with both a tumour and normal sample is treated as one indivisible group. All samples from that patient must remain in the same training, validation or test partition.

## Expression preprocessing

The matrix builder reads the `gene_id` and `tpm_unstranded` columns from each augmented STAR gene-count file, skips GDC summary rows, checks finite non-negative values, aligns genes by stable gene ID, and retains the sample IDs. The baseline applies `log1p(TPM)`.

Within every training fold, the model pipeline performs:
1. Zero-variance feature removal.
2. ANOVA F-test selection of up to 100 genes.
3. Standardization.
4. Class-balanced logistic regression.

Feature selection and scaling are inside the scikit-learn pipeline, so held-out samples do not influence fitted preprocessing.

## Evaluation design

- A patient-grouped 75/25 holdout is generated with a fixed seed.
- Five-fold `StratifiedGroupKFold` cross-validation keeps patient groups intact.
- Metrics: balanced accuracy and ROC-AUC; holdout confusion matrix is also reported.
- The report includes patient IDs in train/test to permit an explicit overlap audit.
- Selected feature coefficients are exploratory model features, not causal biomarkers.

## Interpretation

A strong tumour-versus-normal score on this small, paired sample set would only show that the expression pipeline can recover known specimen labels under the stated split. It would not establish a clinically useful classifier or validate tumour progression, immune escape, spatial oxygen fields, therapy response or patient outcome prediction.

For the larger ONCO-DYNAMICS X research objective, subsequent analyses must use clearly defined clinical endpoints and independent cohorts. Survival analysis requires time-to-event and censoring metadata; treatment-response prediction requires appropriate treatment exposure and response labels. Spatial claims require spatially resolved, registered tissue observations. Bulk expression alone cannot validate cell positions, cell contacts, oxygen gradients or drug penetration.

## Reproducibility

The workflow runs from the repository commit and records the query constraints, selected metadata, manifest, model configuration and metrics. The raw files and expression matrix are transient runner files and are not uploaded as artifacts. To reproduce the analysis, rerun the workflow at the recorded commit; GDC file UUIDs and MD5 checksums in the selected manifest identify the exact source files.

## Data-access boundary

The workflow selects open-access GDC files only. Controlled-access data are not requested. Do not place restricted or identifiable participant-level information in this repository.
