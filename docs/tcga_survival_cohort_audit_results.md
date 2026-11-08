# TCGA-BRCA PanCancer Atlas survival-cohort audit results

## Audit provenance

- Repository: https://github.com/culprit-k/ONCO-DYNAMICS-X
- Workflow: https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/42453906243
- Workflow conclusion: success
- cBioPortal study: brca_tcga_pan_can_atlas_2018
- Study: Breast Invasive Carcinoma (TCGA, PanCancer Atlas)
- Study citation in cBioPortal metadata: TCGA, Cell 2018
- Reference genome metadata: hg19

The workflow queried the public cBioPortal API and uploaded aggregate metadata only. No patient IDs, sample IDs or individual clinical values are stored in the artifact or this repository document.

## Cohort and sample linkage

| Audit item | Result |
|---|---:|
| Unique patients | 1,084 |
| Unique samples | 1,084 |
| Samples linked to a known patient | 1,084 |
| Orphan samples | 0 |
| Clinical attributes | 60 |
| mRNA RNA-seq sample-list members | 1,082 |
| Expression samples linked to a patient | 1,082 |
| Expression patients with complete OS fields | 1,082 |
| Expression patients with positive OS follow-up | 1,069 |

The expression profile is mRNA Expression, RSEM, batch-normalized from Illumina HiSeq RNASeqV2. The profile is distinct from the STAR-count/TPM files used by the earlier GDC paired tumour/normal benchmark, so these matrices must not be treated as interchangeable without explicit harmonization.

## Overall survival

The patient-level fields OS_MONTHS and OS_STATUS are both present for all 1,084 patients.

| OS audit item | Result |
|---|---:|
| Complete time/status pairs | 1,084 |
| Deceased (event = 1) | 151 |
| Living (censored = 0) | 933 |
| Positive follow-up times | 1,071 |
| Non-positive follow-up times | 13 |
| Positive follow-up range | 0.0329–282.901 months |
| Conflicting duplicate time/status records | 0 |

All 13 non-positive follow-up times belong to patients recorded as living/censored; all 151 recorded deaths have positive follow-up times. A conservative model-ready OS subset therefore contains 1,071 patients, including all 151 deaths and 920 living/censored patients.

Within the 1,082-sample mRNA subset:
- 1,082 patients have complete OS time/status fields.
- 1,069 have positive OS follow-up.
- Those 1,069 include 151 deaths and 918 living/censored patients.

The cBioPortal field convention describes OS_MONTHS as months from initial diagnosis and OS_STATUS as 1 for death/event and 0 for living/no event. The cBioPortal FAQ cautions that study-specific endpoint definitions may vary, so the TCGA source documentation must still be checked before the endpoint is locked.

## Other time-to-event fields

| Endpoint | Complete pairs | Event/status categories across complete pairs | Non-positive times |
|---|---:|---|---:|
| Disease-specific survival (DSS) | 1,063 | 83 dead with tumour; 980 alive or dead tumour-free | 13 |
| Disease-free survival (DFS) | 941 | 84 recurred/progressed; 857 disease-free | 11 |
| Progression-free survival (PFS) | 1,082 | 145 progression; 937 censored | 13 |
| Relapse-free survival (RFS) | Not available | No RFS field pair returned | — |

All non-positive times in these endpoint pairs are associated with the no-event/censored status category. The source-specific definitions for DSS, DFS and PFS still require review; do not treat these endpoints as interchangeable.

## Candidate baseline covariates

The audit found the following patient-level fields:

| Field | Patients with a value |
|---|---:|
| Age at diagnosis (AGE) | 1,084 |
| AJCC pathologic tumour stage | 1,079 |
| Pathologic T stage | 1,084 |
| Pathologic N stage | 1,084 |
| Pathologic M stage | 1,084 |
| Winter hypoxia score | 1,079 |

Complete-case overlap:
- Age + AJCC stage: 1,079 patients.
- Age + AJCC stage + Winter hypoxia score: 1,074 patients.
- Within the mRNA subset with positive OS follow-up: 1,064 have age + stage, and 1,061 have age + stage + hypoxia score.

The 1,061-patient group is a candidate for a clinical-plus-hypoxia baseline after field definitions, scale, and timing are verified. Treatment-related variables such as neoadjuvant therapy and radiation must not be used as baseline predictors without careful temporal alignment; they may occur after diagnosis and introduce treatment-selection or immortal-time bias.

## External cohort compatibility: METABRIC

The separately audited METABRIC study has:
- 1,980 mRNA expression samples linked to patients.
- 1,980 complete OS time/status pairs within that expression subset.
- 1,979 positive OS follow-up times.
- 1,143 deaths and 836 living/censored patients among those 1,979 positive-time cases.

METABRIC expression is Illumina HT-12 v3 microarray log2 intensity, whereas this TCGA study uses batch-normalized RNA-seq RSEM. A molecular signature cannot be transported by simply concatenating or directly reusing raw expression values. A frozen TCGA model will require a defensible common-gene mapping and prespecified cross-platform normalization, with cohort/batch diagnostics.

A clinical-only age/stage baseline may be evaluated separately, after stage categories and timing are harmonized. That would be external validation of a clinical model, not validation of the molecular digital twin.

## Current status

Established:
- A public TCGA-BRCA PanCancer Atlas cohort with 1,084 linked patients/samples.
- A 1,082-sample RNA-seq expression subset with complete OS fields.
- 1,069 expression patients with positive OS follow-up and all 151 deaths retained.
- Candidate age/stage/hypoxia covariate overlap quantified.
- A separate METABRIC expression/OS cohort with 1,979 positive-time cases.

Not yet done:
- No Cox model has been fitted.
- No TCGA-to-METABRIC external validation has been performed.
- No gene-level expression harmonization has been implemented or validated.
- No clinical utility or treatment-response claim is supported.
- The digital twin's spatial PDE/ABM mechanisms are not validated by this cohort audit.

## Next implementation gate

Build an endpoint-locked survival dataset assembler and a transparent Cox proportional-hazards baseline:
1. Map OS_STATUS to event (1 = deceased) and censoring (0 = living), and use OS_MONTHS in months.
2. Quarantine the 13 non-positive follow-up cases with an auditable exclusion count; do not silently change their time.
3. Build the TCGA development cohort with patient-level grouping and predeclared age/stage covariates; evaluate the hypoxia score as a separately specified extension.
4. Freeze all preprocessing and coefficients before any METABRIC evaluation.
5. Harmonize stage definitions and test the frozen clinical model on METABRIC first. Only then plan molecular cross-platform validation.

## References

- TCGA PanCancer Atlas study metadata: https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018
- cBioPortal FAQ on OS_STATUS / OS_MONTHS: https://docs.cbioportal.org/user-guide/faq/
- cBioPortal clinical field conventions: https://github.com/cbioportal/cbioportal/blob/master/docs/File-Formats.md
- NCI Genomic Data Commons PanCancer Atlas: https://gdc.cancer.gov/about-data/publications/pancanatlas
