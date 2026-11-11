# Public breast-cancer molecular marker availability audit

**Workflow:** `.github/workflows/molecular-audit.yml`  
**Script:** `scripts/run_molecular_marker_audit.py`  
**Implementation:** `src/onco_dynamics_x/data_engine/molecular_audit.py`  
**Status:** implementation and CI run pending first complete external API evaluation.

## Aim

Establish whether a small, predeclared descriptive panel has measurements in the selected TCGA-BRCA PanCancer Atlas RNA-seq and METABRIC mRNA expression profiles. This is a data-availability and distribution audit, not a molecular signature or predictive model.

## Cohorts and molecular profiles

- TCGA-BRCA PanCancer Atlas: `brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna`
- METABRIC: `brca_metabric_mrna`

The runner checks the study, molecular-profile and sample-list metadata before requesting expression data. It retrieves each marker from the matching study profile and sample list. Gene identifiers are checked against the cBioPortal gene table before requests are accepted.

The cBioPortal public API documents profile-level molecular data retrieval by molecular profile, sample list and Entrez gene identifier, and exposes the response as numeric molecular data. API documentation is marked beta and may change; all requests are validated and a failed or mismatched profile stops the audit.

## Fixed descriptive panel

| Marker | Entrez ID | Context label |
|---|---:|---|
| ESR1 | 2099 | Hormone receptor signalling |
| PGR | 5241 | Hormone receptor signalling |
| ERBB2 | 2064 | Receptor tyrosine kinase signalling |
| MKI67 | 4288 | Proliferation |
| TP53 | 7157 | Tumour suppressor |
| BRCA1 | 672 | DNA repair |
| BRCA2 | 675 | DNA repair |
| EGFR | 1956 | Receptor tyrosine kinase signalling |
| HIF1A | 3091 | Hypoxia response |
| VEGFA | 7422 | Angiogenesis |
| CA9 | 768 | Hypoxia response |
| SLC2A1 | 6513 | Glucose transport |
| STAT1 | 6772 | Interferon signalling |
| CD274 | 29126 | Immune checkpoint |
| CXCL9 | 4283 | Immune chemokine |
| CD8A | 925 | T-cell marker |

These are descriptive pathway/context labels for data availability only. The panel is not a validated signature and has not been selected or tuned using survival, response, or other clinical outcomes.

## Per-marker output

For each gene and profile, the report contains:
- returned record count;
- unique samples with an expression record;
- fraction of the predefined sample list covered;
- missing, non-numeric/non-finite, duplicate and identifier-free record counts;
- finite value count;
- within-profile min, p10, p25, median, p75, p90, max, mean and standard deviation.

Cross-profile overlap indicates the number of panel genes with finite values in every queried profile. The audit does not compare the numeric distributions directly between assays.

## Important cross-platform boundary

The audited TCGA PanCancer Atlas expression profile is RNA-seq RSEM-derived; METABRIC expression is Illumina HT-12 v3 microarray log2 intensity. These measurements have different assay technologies, preprocessing and scales. They must not be concatenated, z-scored together without a prespecified plan, or treated as the same quantity. This audit is a feasibility check only.

A later molecular-model pipeline would still need:
1. review of the original cohort assay/probe annotation and version;
2. versioned Entrez/HUGO/probe mapping, including duplicate probe/gene rules;
3. a predeclared gene set and biological question;
4. a defensible training-only cross-platform transformation;
5. missing-gene and batch diagnostics;
6. an outcome-locked TCGA development design and frozen independent METABRIC evaluation;
7. discrimination, calibration and uncertainty assessment.

No external outcome labels are queried in this workflow.

## Data governance

Only public cBioPortal expression profiles are queried. Raw molecular records remain in memory and are reduced to aggregate summaries. Patient/sample IDs and raw expression values are not written to files or uploaded as workflow artifacts. The workflow publishes only the aggregate JSON report, retained for 30 days.

## References

- cBioPortal API and clients: https://github.com/cBioPortal/cbioportal/blob/master/docs/web-API-and-Clients.md
- cBioPortal current OpenAPI specification: https://www.cbioportal.org/api/v3/api-docs
- cBioPortal clinical and molecular data conventions: https://docs.cbioportal.org/file-formats/
- cBioPortalData molecular-data workflow: https://waldronlab.io/cBioPortalData/reference/cBioPortal.html

## Interpretation boundary

The audit can establish marker data availability, coverage, missingness and within-profile distribution summaries. It cannot establish:
- a validated prognostic or diagnostic signature;
- treatment-response prediction or treatment benefit;
- causal roles of any marker;
- hypoxia or immune-cell spatial distributions;
- tumour dynamics or drug penetration;
- biological or clinical validity of the ONCO-DYNAMICS-X digital twin.
