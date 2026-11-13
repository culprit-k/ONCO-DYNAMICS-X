# I-SPY 2 252-participant DICOM completion workflow

## Verified metadata-level status

The Colab audit reported:

- 384 unique participants in the feature workbook.
- 719 MR PatientIDs in the TCIA ISPY2 collection.
- 697 successfully retrieved ClinicalTrialSubjectID DICOM values.
- 273 participants matched to the feature workbook.
- Matched-cohort visit availability: T0 270, T1 264, T2 266, T3 268.
- 252 matched participants with all four T0–T3 visits.
- 1,008 selected original-DCE source series, one per matched participant and visit.
- 803 unilateral-cropped and 205 bilateral original-DCE descriptions.

These are metadata audit results only. They do not establish that images have been downloaded, pixel data inspected, maps registered, FTV recomputed or biological parameters calibrated. The 252-person set remains a candidate longitudinal imaging cohort.

## Workflow

### 1. Expand each selected study before downloading images

The original-DCE source series may not contain the entire DCE observation package. PE maps, SER maps, the bit-encoded analysis mask and anatomical reference series can be separate series within the same examination.

Run the new scripts/expand_ispy2_dce_study_series.py utility against the private 1,008-row review CSV. It queries TCIA v4 one PatientID at a time, filters to the selected StudyInstanceUIDs, caches each patient response and creates a private candidate table. Candidate rows always default to include=0. This script does not download images or choose series.

Review candidate role counts, metadata and image content. Set include=1 only for series actually reviewed, and replace the default review note with a real rationale. Use the existing reviewed-manifest exporter after review. Series descriptions alone do not prove correct image content, geometry or anatomical suitability.

### 2. Download only the reviewed manifest

The new scripts/download_tcia_manifest_with_retriever.py utility downloads the latest official TCIA Data Retriever Linux x86_64 release, verifies the SHA-256 digest published with the release asset, caches the executable and invokes it with server-friendly and skip-existing options. This allows interrupted downloads to be resumed.

The helper does not silently accept the TCIA data-use policy. Read the current policy and accept it yourself. Only after personal acceptance, and only if the non-interactive run requires it, pass the explicit accept-data-policy option. Do not store DICOMs, manifest files, logs or participant-level data in the public repository.

### 3. Inventory and build the patient timeline

The existing repository tools include:
- scripts/inventory_ispy2_dicoms.py for recursive DICOM header inventory without reading pixel arrays;
- scripts/build_ispy2_visit_timeline.py for DICOM visit and date-provenance checks;
- scripts/register_ispy2_dicom_series.py for first-pass rigid anatomical registration;
- scripts/resample_ispy2_dce_map.py for map resampling with explicit frame-of-reference checks;
- extended-Tofts/SPGR components and an FTV measurement-side operator.

The timeline tool blocks ambiguous repeat studies, missing or inconsistent date-modification flags, non-monotone order and inadequate time precision. TCIA-shifted dates are pseudodates: use only defensible within-participant elapsed intervals, never infer real calendar dates.

### 4. Required image-level gates

Before quantitative biological modelling:

1. Verify downloaded SeriesInstanceUID coverage against the reviewed manifest; inspect missing, failed and partial series.
2. Review all four studies per participant, including unilateral/bilateral coverage and optional retest ambiguity.
3. Obtain PEearly, PElate, SER, analysis-mask and suitable anatomical-reference series.
4. Verify rescale slope/intercept, voxel spacing, orientation, frame of reference, sequence parameters and image dimensions.
5. Obtain and verify the official FTV mask bit assignments. Never treat the raw encoded mask as binary using mask > 0.
6. Recompute FTV only with verified mask semantics and thresholds; compare with source DICOM FTV result values as a QC check.
7. Review visual overlays and landmarks after registration. Rigid registration is not proof of breast-deformation correction.
8. Obtain actual administered treatment events and a documented mapping to the TCIA modified date axis. Treatment-arm labels are not actual dose histories.
9. Validate the latent PDE-to-MRI observation operator and parameter identifiability before estimating biological parameters.
10. Split train/validation/test by participant so that all scans from one participant remain in one partition.

No proliferation, cell-death, drug-effect, immune-killing or diffusion coefficient should be described as biologically calibrated until measurement, timing, exposure and identifiability gates pass.

## Privacy

Keep participant-level review tables, IDs, DICOM UIDs, shifted dates, metadata caches, DICOMs, local paths, maps, masks, registration transforms and per-participant QC in approved private storage. Do not commit them to this public repository. Only privacy-reviewed aggregate summaries are eligible for publication.

## Colab notebook

A guided notebook is provided separately as the ONCO-DYNAMICS-X_ISPY2_252_DICOM_pipeline notebook. It mounts Google Drive before creating persistent directories, revalidates the 252-participant and 1,008-series manifest, expands sibling series, gates series selection on human review, and provides download, inventory and timeline steps.
