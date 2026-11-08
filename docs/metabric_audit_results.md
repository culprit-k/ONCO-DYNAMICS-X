# METABRIC cohort audit results

## Audit provenance

- Repository: https://github.com/culprit-k/ONCO-DYNAMICS-X
- Workflow: https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/42453381822
- Workflow result: success
- cBioPortal study: brca_metabric
- Study release presented by cBioPortal: Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016)
- Reference genome metadata: hg19
- Study citations in cBioPortal metadata: Curtis et al. (Nature 2012), Pereira et al. (Nature Communications 2016), and Rueda et al. (Nature 2019).

The audit is reproducible through the METABRIC independent cohort audit GitHub Actions workflow. The workflow artifact contains aggregate metadata only. It does not contain patient IDs, sample IDs, or individual clinical values.

## Cohort and linkage

| Audit item | Result |
|---|---:|
| Unique patients | 2,509 |
| Unique samples | 2,509 |
| Samples linked to a known patient | 2,509 |
| Orphan samples | 0 |
| Clinical attributes | 36 |
| mRNA expression sample-list members | 1,980 |
| Expression samples linked to a patient | 1,980 |
| Expression patients with complete OS fields | 1,980 |
| Expression patients with complete RFS fields | 1,979 |

The cBioPortal study description identifies 2,509 primary breast tumours and 548 matched normal samples. The mRNA sample list contains 1,980 samples. The patient-to-sample mapping was complete for all 2,509 samples returned by the study endpoint.

## Molecular profile inventory

Five molecular profiles were returned:

| Profile | Assay / data type |
|---|---|
| mRNA expression | Illumina HT-12 v3 microarray, log2 intensity |
| mRNA expression z-scores | Log2-transformed microarray z-scores relative to all samples |
| Copy-number alterations | DNAcopy-derived discrete CNA calls |
| Promoter methylation | RRBS beta-values |
| Mutations | MAF from targeted sequencing of 173 genes |

**Critical cross-cohort issue:** METABRIC expression is microarray log2 intensity. The current TCGA GDC pipeline uses STAR-count-derived TPM. These are not interchangeable scales or platforms. Do not concatenate them or apply the TCGA tumour/normal classifier directly to METABRIC. Any molecular external validation needs a prespecified, defensible gene mapping and platform-harmonization strategy, followed by cohort/batch diagnostics.

## Overall-survival endpoint audit

The patient-level clinical data contain:

- OS_MONTHS: 1,981 non-empty values.
- OS_STATUS: 1,981 non-empty values.
- Patients with both fields: 1,981.
- Conflicting duplicate values: 0 for both fields.
- Numeric follow-up values: 1,981.
- Positive follow-up values: 1,980.
- Zero or negative follow-up values: 1.
- Positive follow-up range: 0.1–355.2 months.
- Status values: 1,144 deceased and 837 living.

The cBioPortal clinical-data convention describes OS_MONTHS as months since initial diagnosis and OS_STATUS as 1 for an event/death and 0 for no event/living. The portal documentation also cautions that endpoint definitions can vary by study. Therefore, the original METABRIC study documentation must still be checked before final endpoint locking. The single non-positive time must be quarantined or resolved using the source record before fitting a model.

The mRNA sample list contains 1,980 samples, and all 1,980 corresponding patients have both OS fields. This establishes sample-level availability/linkage for an OS analysis; it does not establish clinical validity or external model performance.

## Relapse-free endpoint audit

The portal contains:

- RFS_MONTHS: 2,388 non-empty values.
- RFS_STATUS: 2,488 non-empty values.
- Patients with both fields: 2,380.
- Conflicting duplicate values: 0 for both fields.
- Numeric follow-up values among complete pairs: 2,380.
- Positive follow-up values: 2,376.
- Zero or negative follow-up values: 4.
- Positive follow-up range: 0.0333–389.3333 months.

The status categories across all available RFS_STATUS values are 1,002 “Recurred” and 1,486 “Not Recurred”. These counts are for all non-empty status records, not just the 2,380 complete pairs. The current audit does not establish the exact event definition, time origin, censoring convention, or whether death is included as an event in this study-specific RFS endpoint. Keep RFS exploratory until the original data guide/publication resolves those points.

Of the 1,980 patients with mRNA expression data, 1,979 have both RFS fields.

## Treatment fields

Patient-level fields are present for:
- CHEMOTHERAPY: 1,980 records with values.
- HORMONE_THERAPY: 1,980 records with values.
- RADIO_THERAPY: 1,980 records with values.

These are treatment-exposure fields, not treatment-response labels. They do not by themselves establish treatment line, timing, dose, regimen, response criteria, or treatment benefit. No response model is justified by these fields alone.

## What this audit establishes — and what it does not

Established:
- The public study endpoint is reachable.
- The study has 2,509 patients and 2,509 sample records with complete sample-to-patient linkage.
- There is a 1,980-sample mRNA expression subset.
- The expression subset links to 1,980 patients with complete OS fields.
- OS and RFS status/time field pairs exist with the aggregate completeness and data-quality characteristics above.
- Treatment-exposure attributes exist.

Not established:
- A survival model has been fitted.
- A frozen TCGA-trained model has been evaluated on METABRIC.
- Molecular platform harmonization has been solved.
- RFS semantics have been fully verified.
- Any treatment-response or causal treatment-effect model is supported.
- The digital twin's spatial PDE/ABM mechanisms are biologically validated.

## Next gate

Audit a suitable TCGA-BRCA development cohort with patient-linked OS fields and a compatible expression profile. The existing 48-patient paired TCGA tumour/normal benchmark is not a suitable survival-model training cohort. Once a larger TCGA cohort and its endpoint linkage are verified, lock a TCGA-only development pipeline, harmonize the molecular feature space without using METABRIC outcomes, and evaluate the frozen model on the 1,980-patient METABRIC expression/OS overlap.

## References

- cBioPortal FAQ: https://docs.cbioportal.org/user-guide/faq/
- cBioPortal clinical data file format and OS field conventions: https://github.com/cbioportal/cbioportal/blob/master/docs/File-Formats.md
- METABRIC study: https://www.cbioportal.org/study/summary?id=brca_metabric
