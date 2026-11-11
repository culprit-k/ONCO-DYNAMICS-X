# Phase 3 — fixed-panel molecular model development

**Workflow:** `.github/workflows/pam50-molecular-validation.yml`  
**Runner:** `scripts/run_pam50_molecular_validation.py`  
**Feature implementation:** `src/onco_dynamics_x/data_engine/molecular_model.py`  
**Status:** code and tests committed; report is produced by the external-cohort workflow.

## Research question

Does a compact, predeclared breast-cancer expression panel contribute measurable prognostic information beyond age at diagnosis when a model developed in TCGA-BRCA is transported to METABRIC?

This is an exploratory retrospective prognostic-association experiment. It does not ask whether a treatment works, does not estimate an individual patient's risk for clinical use, and does not validate the spatial PDE/agent-based digital twin.

## Fixed feature inventory

The input panel is the 50-gene intrinsic-subtype research panel exposed by the `genefu` PAM50 model mapping:

```
ACTR3B, ANLN, BAG1, BCL2, BIRC5, BLVRA, CCNB1, CCNE1,
CDC20, CDC6, CDCA1, CDH3, CENPF, CEP55, CXXC5, EGFR,
ERBB2, ESR1, EXO1, FGFR4, FOXA1, FOXC1, GPR160, GRB7,
KIF2C, KNTC2, KRT14, KRT17, KRT5, MAPT, MDM2, MELK,
MIA, MKI67, MLPH, MMP11, MYBL2, MYC, NAT1, ORC6L,
PGR, PHGDH, PTTG1, RRM2, SFRP1, SLC39A6, TMEM45B,
TYMS, UBE2C, UBE2T
```

The panel is fixed before analysis. Historical symbols are resolved through the cBioPortal gene catalogue using current-symbol aliases (for example, CDCA1/NUF2, KNTC2/NDC80 and ORC6L/ORC6); both the requested and resolved symbols are retained in aggregate output.

**Important distinction:** this project does not import or reproduce the published PAM50 centroids, the commercial Prosigna assay, its proprietary Risk-of-Recurrence algorithm, clinical cut-points or intended-use interpretation. The panel is used only as a published, fixed gene inventory for transparent research feature construction. A five-PC model built from it is not a PAM50 intrinsic-subtype classifier.

## Cohorts and data

- Development: TCGA-BRCA PanCancer Atlas RNA-seq mRNA profile `brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna`.
- External evaluation: METABRIC mRNA profile `brca_metabric_mrna`.
- Samples are restricted to primary tumour where sample-type metadata support that filtering. If sample-type metadata are unavailable from the public API, that limitation is explicitly reported and the profile sample list is retained.
- One expression specimen per patient is included in the survival analysis. Repeated specimens are counted and excluded from the primary patient-level analysis.
- Overall survival uses `OS_MONTHS` and `OS_STATUS`; 1:DECEASED is an event and 0:LIVING is censored under the cBioPortal convention. Non-positive follow-up is excluded; no artificial time is substituted.
- Stage is not included because the TCGA patient-level AJCC pathologic stage and METABRIC sample-level generic tumour-stage field have not been established as semantically equivalent.

The exact time origin and endpoint adjudication remain study-specific and must be checked against the original cohort publications/data guides before any clinical interpretation.

## Outcome-blind expression processing

### 1. Gene resolution and completeness audit

Each panel symbol is resolved against the cBioPortal gene catalogue. Expression is fetched in bounded batches using sample IDs and Entrez identifiers. The runner records aggregate counts for missing, duplicate, non-finite and unrecognised sample-gene records. It never writes sample or patient identifiers or raw expression values.

A gene is retained only when at least 95% of samples have finite values in **both** cohorts. At least 45 shared genes are required for the analysis to proceed. The actual retained gene list and coverage by gene are reported.

### 2. Cohort-wise rank-normalisation

Within each cohort independently, each gene is transformed using its within-cohort mid-ranks:

```
p_ig = (rank_ig - 0.5) / n_g
z_ig = inverse_normal_CDF(p_ig)
```

Tied values receive average ranks. Missing values, if coverage remains above the prespecified threshold, are assigned zero on this rank-normal scale (the cohort median). No endpoint, subtype annotation or treatment outcome enters this step.

This is a **transductive, cohort-level** harmonisation: all unlabeled expression values in the METABRIC evaluation cohort contribute to its per-gene ranking. It can reduce between-platform marginal scale differences but cannot guarantee removal of platform bias, probe effects, batch effects or differing sample composition. It is not a deployable single-patient normaliser. A future prospective pipeline would need a separately specified reference distribution or an assay-specific frozen calibration.

### 3. TCGA-fitted principal components

The rank-normalised TCGA matrix is used to fit five-component PCA. The resulting component loadings are frozen and applied unchanged to the rank-normalised METABRIC matrix. The model report includes explained variance and the highest positive/negative loadings for each component.

The PCA is unsupervised with respect to survival and subtype labels. It does not infer biological pathways from gene labels automatically; its components are statistical axes whose biological meaning requires separate validation.

## Survival models and comparison

Two penalised Cox models are fitted in TCGA using the same patient set:

1. **Reference:** age at diagnosis only.
2. **Molecular model:** age plus the five TCGA-fitted expression principal components.

Age centring/scaling is learned in TCGA and carried forward unchanged. The molecular model uses an L2 penalizer of 0.5. The external METABRIC outcomes are held out from feature construction and model fitting.

Primary external measures:
- Harrell's C-index with patient-bootstrap confidence intervals;
- paired bootstrap difference in C-index versus the age-only comparator;
- cumulative/dynamic time-dependent AUC at 12, 36 and 60 months;
- IPCW Brier scores at the same horizons;
- descriptive IPCW calibration summaries;
- proportional-hazards diagnostic.

The IPCW procedure estimates a marginal censoring Kaplan-Meier curve within METABRIC for evaluation-weight correction. This assumes independent censoring in the external cohort. Covariate-dependent censoring is not modelled. Metrics that cannot be estimated because of insufficient cases/controls or unstable censoring survival are reported as not estimable rather than fabricated.

No threshold is selected using METABRIC survival. A performance increase over age-only is not automatically clinically meaningful.

## Subtype annotation handling

The runner inventories available TCGA/METABRIC clinical attributes that appear to contain PAM50 or molecular-subtype annotations. It records aggregate value counts only. These annotations are not used in gene selection, normalization, PCA, model fitting or hyperparameter tuning. If sufficiently reliable annotations are available in both cohorts, a separate descriptive concordance study can be prespecified later. This run does not claim intrinsic subtype classification.

## Data governance

Only aggregate counts, distributions, model coefficients, metrics and PCA loadings are written to the report. No patient IDs, sample IDs, raw expression values, or row-level clinical records are written or uploaded. Controlled/restricted source data must remain in approved storage.

## Reproduction

```bash
pip install -e ".[dev,science,survival]"
pytest
ruff check src tests scripts/run_pam50_molecular_validation.py
python scripts/run_pam50_molecular_validation.py \
  --output results/molecular_validation/pam50_tcga_metabric_survival.json \
  --components 5 \
  --bootstrap 250 \
  --seed 42 \
  --cox-penalizer 0.5
```

The GitHub Actions workflow runs the complete test and lint suite, executes the public cohort evaluation and publishes only the aggregate JSON report.

## Interpretation limits

This experiment is retrospective and exploratory. It does not establish:
- a validated PAM50 or Prosigna clinical test;
- a clinically deployable single-patient predictor;
- predictive response to endocrine therapy, chemotherapy, HER2-targeted therapy or immunotherapy;
- a causal role for any panel gene;
- a spatial immune-infiltration or tumour-microenvironment measurement;
- validation of PDE oxygen gradients, cell contacts, drug penetration or digital-twin trajectories;
- clinical utility or improved patient outcomes.

Any molecular association must be interpreted alongside the external cohort's case mix, censoring, data provenance, platform differences and model calibration.

## References

- Parker JS et al. *Supervised Risk Predictor of Breast Cancer Based on Intrinsic Subtypes*. Journal of Clinical Oncology. 2009;27:1160–1167. https://doi.org/10.1200/JCO.2008.18.1370
- `genefu` PAM50 gene mapping and centroid metadata (gene inventory source; centroids are not copied into this project): https://github.com/bhklab/genefu/blob/master/inst/extdata/pam50_model.csv
- `genefu` PAM50 model documentation: https://rdrr.io/bioc/genefu/man/pam50.html
- cBioPortal API and clients: https://github.com/cBioPortal/cbioportal/blob/master/docs/web-API-and-Clients.md
- cBioPortal OpenAPI specification: https://www.cbioportal.org/api/v3/api-docs
- cBioPortal clinical and molecular data conventions: https://docs.cbioportal.org/file-formats/
