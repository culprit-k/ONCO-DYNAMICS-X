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


## Completed audit — 11 November 2026

**Workflow:** [run 42982626324](https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/42982626324)  
**Evaluated commit:** `1f251ed54031130909ae9dd62808b0306055c2f0`  
**Committed aggregate report:** `results/molecular_audit/marker_panel_overlap_summary.json`

| Cohort | Profile | Sample count | Markers with finite values | Marker-sample values returned | Missing / non-finite | Duplicate sample records |
|---|---|---:|---:|---:|---:|---:|
| TCGA-BRCA PanCancer Atlas | RSEM batch-normalized RNA-seq | 1,082 | 16 / 16 | 17,312 / 17,312 | 0 | 0 |
| METABRIC | Illumina HT-12 v3 microarray log2 intensity | 1,980 | 16 / 16 | 31,680 / 31,680 | 0 | 0 |

Each marker had 100% sample-list coverage in both profiles. All Entrez IDs resolved to the expected HUGO symbols. No patient or sample identifiers or raw expression values were written to the report.

### Selected within-profile summaries

| Gene | TCGA median (p25–p75) | METABRIC median (p25–p75) |
|---|---:|---:|
| ESR1 | 7,065.11 (1,251.32–15,618.05) | 10.2235 (7.9225–11.2588) |
| PGR | 726.782 (58.526–3,870.138) | 5.8665 (5.4083–6.8968) |
| ERBB2 | 7,025.295 (4,476.123–10,792.15) | 10.5260 (9.9584–11.1618) |
| MKI67 | 1,699.535 (952.002–2,981.413) | 5.8125 (5.6159–6.0551) |
| HIF1A | 2,159.915 (1,296.903–3,446.438) | 6.6309 (6.2069–7.0939) |
| VEGFA | 1,753.415 (1,212.743–2,704.330) | 7.0711 (6.6975–7.5543) |
| CA9 | 3.9656 (0.5407–47.4461) | 5.5828 (5.4417–5.9143) |
| CD274 | 19.7097 (10.8313–37.4860) | 5.4751 (5.3800–5.5782) |
| CXCL9 | 690.9985 (234.095–2,464.820) | 8.2036 (6.9615–9.6484) |
| CD8A | 153.8000 (73.7001–325.5670) | 7.1708 (6.4483–7.9641) |

The difference in numeric scale is conspicuous and expected from the distinct measurement platforms and processing. The table is for within-cohort descriptive reporting only. **It does not imply that values from TCGA and METABRIC are on a common scale or directly comparable.**

### Decision from this audit

The selected marker panel is sufficiently represented in both cohorts to justify a carefully designed next-stage harmonization feasibility study. This does not establish that an expression signature is transportable. Before fitting a molecular predictor, the next required work is probe-to-gene mapping for the METABRIC microarray, duplicate-probe summarization, versioned gene identifier mapping, a training-only transformation, and a frozen external evaluation design. No METABRIC outcome should be used to choose genes, transformations, cut-offs or model parameters.

**Current status:** expression availability verified for this small marker panel; molecular harmonization and molecular model training remain pending.
