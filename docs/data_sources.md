# Public cohort sources and GDC workflow

## NCI Genomic Data Commons

The repository contains a small GDC API client for public/open file metadata. The NCI GDC API supports filtered search and retrieval through endpoints including files, cases and projects. Official references:
- [GDC API Search and Retrieval](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/)
- [GDC API Getting Started](https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/)
- [GDC API overview](https://gdc.cancer.gov/developers/gdc-application-programming-interface-api)
- [GDC Data Access Processes and Tools](https://gdc.cancer.gov/access-data/data-access-processes-and-tools)

The default example is TCGA-BRCA and Gene Expression Quantification. No workflow type is assumed by the client unless the caller supplies one.

### Step 1: query public metadata

```bash
onco-dynamics gdc-query \
  --project-id TCGA-BRCA \
  --data-type "Gene Expression Quantification" \
  --workflow-type "STAR - Counts" \
  --max-results 100
```

This command uses an explicit `open` access filter and writes local metadata and sample-manifest CSV files. The metadata records include GDC file identifiers, filename, size, checksum, workflow, case and sample identifiers. The manifest is audited before analysis. It queries metadata only and does not download files.

### Step 2: review the metadata and sample linkage

Review `data/metadata/gdc_open_files.csv` and `data/metadata/gdc_sample_manifest.csv`. Verify that file/sample relationships, sample types, assay type, project, analysis workflow and intended endpoint are suitable for the intended research question. The manifest uses source submitter identifiers for auditability; these are not patient medical records but should still be handled under data-use and institutional policy. Do not publish any row-level data unless access terms explicitly allow it.

### Step 3: download one selected open file

After review, choose a specific file UUID from the metadata table:

```bash
onco-dynamics gdc-download-open \
  --metadata-csv data/metadata/gdc_open_files.csv \
  --file-id YOUR-GDC-FILE-UUID \
  --destination data/raw/gdc
```

The command downloads only the selected metadata row. The client refuses files not marked open, validates the GDC UUID and plain filename, then verifies the exact expected byte count and MD5 checksum before finalizing the download. It does not accept GDC authentication tokens and must not be used to access controlled data. For large-volume downloads, follow the GDC Data Transfer Tool and official access procedures.

### Other potential cohorts

- [METABRIC in cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric): a potential independent breast-cancer cohort. Confirm current study metadata, endpoint definitions, molecular platforms and permitted access before analysis.
- [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/): expression and other assay series. Each series needs a specific accession audit and sample-annotation review.
- [Proteomic Data Commons](https://proteomic.datacommons.cancer.gov/): proteomics source. Confirm cohort and sample linkage before joining with another modality.

No two cohorts or modalities should be treated as paired merely because they concern breast cancer. Preserve patient/specimen/sample linkage from source metadata. Separate prognostic survival modelling from treatment-response or treatment-benefit estimation. Bulk RNA expression does not provide spatial geometry, oxygen gradients or cell-cell contact ground truth.

## Data governance

- Do not upload controlled, restricted or identifiable participant-level data to this public repository.
- Keep raw and intermediate data in approved local or institutional storage.
- Record dataset version, download date, source accession, checksums, inclusion criteria and preprocessing version.
- Respect the relevant GDC, dbGaP, institutional and data-use requirements.
