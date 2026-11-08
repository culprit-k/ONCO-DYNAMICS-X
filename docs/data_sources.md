# Public cohort sources and GDC workflow

## NCI Genomic Data Commons

The repository includes a small GDC API client for public/open file metadata. The NCI GDC API supports filtered search and retrieval through endpoints including files, cases and projects. Its search guide documents POST requests with JSON filters, fields, pagination and response formats:
- [GDC API Search and Retrieval](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/)
- [GDC API Getting Started](https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/)
- [GDC Data Portal](https://portal.gdc.cancer.gov/)

The default example is TCGA-BRCA and Gene Expression Quantification, but the API client does not assume a particular workflow unless one is provided.

### Metadata-only query

```bash
onco-dynamics gdc-query \
  --project-id TCGA-BRCA \
  --data-type "Gene Expression Quantification" \
  --workflow-type "STAR - Counts" \
  --max-results 100
```

This command:
1. queries file metadata using an explicit `open` access filter;
2. saves a local metadata CSV with file IDs, filenames, size, checksums, workflow and case/sample identifiers;
3. creates a sample-level manifest;
4. audits that manifest using the repository's required-column and identifier checks.

Generated CSVs are ignored by Git by default. Review them locally before sharing. The command queries metadata only; it does not download hundreds of files by default.

### Download policy

The API client exposes a separate method for an individual metadata record. It refuses non-open files, requires a valid GDC UUID, checks the expected byte count and verifies the MD5 checksum before moving the temporary file into place. This is not a bulk downloader. Do not try to bypass controlled-access requirements. Follow the GDC/NIH data-use and institutional rules for any controlled data; never commit controlled or identifiable data to this repository.

## Other potential cohorts

- [METABRIC in cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric): a potential independent breast-cancer cohort. Confirm current study metadata, endpoint definitions, molecular platforms and permitted access before analysis.
- [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/): expression and other assay series. Each series needs a specific accession audit and sample-annotation review.
- [Proteomic Data Commons](https://proteomic.datacommons.cancer.gov/): proteomics source. Confirm cohort and sample linkage before joining with another modality.

No two cohorts or modalities should be treated as paired merely because they concern breast cancer. Preserve patient and specimen linkage from source metadata. Separate prognostic survival modelling from treatment-response or treatment-benefit estimation.
