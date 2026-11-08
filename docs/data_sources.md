# Public cohort sources and GDC workflow

## NCI Genomic Data Commons

The repository contains a standard-library GDC API client for querying public/open file metadata and downloading a specifically selected open file. The NCI GDC API supports filtered search and retrieval through endpoints including files, cases and projects. Official references:
- [GDC API Search and Retrieval](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/)
- [GDC API Getting Started](https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/)
- [GDC API overview](https://gdc.cancer.gov/developers/gdc-application-programming-interface-api)
- [GDC Data Access Processes and Tools](https://gdc.cancer.gov/access-data/data-access-processes-and-tools)

The example query defaults to project `TCGA-BRCA` and `Gene Expression Quantification`, with an optional analysis workflow filter.

### Step 1: query metadata only

```bash
onco-dynamics gdc-query \
  --project-id TCGA-BRCA \
  --data-type "Gene Expression Quantification" \
  --workflow-type "STAR - Counts" \
  --max-results 100
```

This command restricts results to records marked `open`, saves file metadata and case/sample identifiers, exports a sample-level manifest and audits its required fields. It does not download files.

### Step 2: review sample and file metadata

Review `data/metadata/gdc_open_files.csv` and `data/metadata/gdc_sample_manifest.csv`. Check that the samples are the desired primary tumour/normal/metastatic specimen types, that the assay workflow and expression scale suit the question, and that patient/sample relationships and endpoint definitions are appropriate.

Each assay file must be linked to the correct sample. A patient can contribute multiple samples; keep patient grouping for all train/test splits. The repo provides `group_train_test_split` to keep samples from the same patient together. This is group splitting, not stratified splitting; event prevalence and group counts still need review.

### Step 3: download one selected open file

After reviewing a file row, select its GDC file UUID:

```bash
onco-dynamics gdc-download-open \
  --metadata-csv data/metadata/gdc_open_files.csv \
  --file-id YOUR-GDC-FILE-UUID \
  --destination data/raw/gdc
```

The downloader refuses non-open records, validates the UUID and plain filename, and checks exact file size and MD5 before finalizing. It does not accept credentials or retrieve controlled data.

### Step 4: assemble an expression matrix

After the selected files have all been downloaded into the same directory:

```bash
onco-dynamics gdc-build-expression-matrix \
  --metadata-csv data/metadata/gdc_open_files.csv \
  --data-directory data/raw/gdc \
  --gene-id-column gene_id \
  --expression-column tpm_unstranded \
  --output data/processed/tcga_brca_expression.npz
```

The parser reads the named columns, skips GDC summary rows beginning with `N_`, rejects duplicate gene IDs, rejects non-finite or non-numeric values, requires exactly one sample ID per file and requires an identical gene set across files. It reorders gene columns to match the first file and retains sample/gene identifiers in the compressed NPZ output. Verify actual file headers and the appropriate expression scale for the selected GDC workflow before running this step.

### Other potential cohorts

- [METABRIC in cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric): a potential independent breast-cancer cohort. Confirm study metadata, endpoint definitions, molecular platforms and access terms before analysis.
- [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/): expression and other assay series. Each accession needs a specific data audit and sample-annotation review.
- [Proteomic Data Commons](https://proteomic.datacommons.cancer.gov/): a proteomics source. Confirm cohort and sample linkage before joining with another modality.

Do not join modalities merely because they concern breast cancer. Preserve sample/patient identity from source metadata. Bulk RNA data cannot supply spatial geometry, oxygen gradients or cell-contact ground truth. A survival association is prognostic evidence, not automatically evidence of treatment benefit.

## Data governance
- Do not upload controlled, restricted or identifiable participant-level data to this public repository.
- Keep raw/intermediate data in approved local or institutional storage.
- Record dataset versions, download dates, access class, checksums, inclusion criteria and preprocessing.
- Respect current GDC, dbGaP, institutional and data-use requirements.
