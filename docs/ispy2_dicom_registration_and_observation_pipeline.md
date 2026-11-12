# I-SPY 2 DICOM registration, timing and MRI observation pipeline

**Status: 12 November 2026.** The repository now contains local-only DICOM metadata, longitudinal-date, registration, treatment-event and DCE-to-FTV tooling. No I-SPY 2 DICOM images or TCIA series manifests were present in the local analysis workspace during this implementation. Therefore no participant-specific DICOM timestamps, registration transforms, DCE maps, treatment-exposure alignments or patient-derived PDE parameters have been generated. The scripts are implementation infrastructure, not a completed imaging analysis.

## 1. Available data and unresolved inputs

The linked research workbooks provide 384 participants with four visit-stage MRI feature summaries (T0–T3), and a separate 985-row clinical table with Patient_ID, treatment Arm, HR, HER2, pCR and selected demographic variables. All 384 MRI participants matched one-to-one through the documented join key CLINICAL-TRIAL-SUBJECT-ID = Patient_ID.

The 985-row clinical workbook does **not** contain actual drug administration dates, doses, dose reductions, missed administrations or treatment start/stop timestamps. The treatment-arm label is not a substitute for exposure history.

The DICOM archive is a separate resource. Do not assume its participant count, visit coverage, identifiers, series availability or measurement implementation are identical to the 384-row workbook. Join using the documented six-digit DICOM ClinicalTrialSubjectID tag (0012,0040), after independently validating the match. The trial visit tag is (0012,0050) and its description is (0012,0051).

## 2. Obtain a selected DICOM subset

Use the official [TCIA I-SPY 2 collection page](https://www.cancerimagingarchive.net/collection/ispy2/) and linked ACRIN-6698/I-SPY2 imaging resources to identify the relevant DICOM exams. The archive includes original T2-weighted, DWI and DCE data plus derived maps and masks; the full collection is large. Do not download the entire collection merely to begin this project.

TCIA is migrating its former NBIA workflow to the [TCIA Data Retriever](https://github.com/TCIA/data-retriever). The current Retriever accepts legacy .tcia manifests and newer .csv/.s5cmd manifests. Its CLI supports metadata-only retrieval (--meta) and SeriesInstanceUID CSV input. Use the TCIA portal's patient/series query to constrain the initial selection to the 384-person feature cohort and required visit/series types. Confirm that the portal patient identifier corresponds to the documented trial identifier before exporting a manifest. If the portal's patient-ID field does not expose the Clinical Trial Subject ID, do not assume the values are interchangeable; use the TCIA-provided query/manifest mapping and verify the DICOM header after download.

Suggested acquisition strategy:
1. Query candidate participant IDs and visits through TCIA; preserve the exact exported manifest and query date locally.
2. Review series descriptions and numbers. Select original T2-weighted anatomical series and the DCE derived package needed for PE/SER/FTV, not every reconstruction.
3. Download the selected patient/series manifest using the current TCIA Data Retriever. Keep the manifest, logs and downloaded files in approved private storage.
4. Run the metadata-only inventory script below on downloaded DICOM files. It reads headers with stop_before_pixels=True; it does not load pixel arrays.
5. Reconcile every DICOM subject ID and visit code to the clinical/MRI workbooks before any model fit. Report unmatched participants explicitly.

The public [NBIA v4 API guide](https://wiki.cancerimagingarchive.net/pages/viewpage.action?pageId=70230731) documents public getPatientStudy and getSeries methods. It queries the DICOM PatientID field, which must not be presumed identical to the Clinical Trial Subject ID. Verify the ID mapping rather than sending workbook six-digit values blindly to a PatientID filter.

## 3. DICOM tags and time handling

The metadata parser reads:
- ClinicalTrialSubjectID (0012,0040) — six-digit trial participant identifier;
- ClinicalTrialTimePointID (0012,0050) — raw visit code, with canonical T0/T1/T2/T3 stage extracted while retaining any sub-visit suffix;
- ClinicalTrialTimePointDescription (0012,0051);
- StudyDate/StudyTime (0008,0020)/(0008,0030), SeriesDate/SeriesTime (0008,0021)/(0008,0031), AcquisitionDate/AcquisitionTime (0008,0022)/(0008,0032) and AcquisitionDateTime (0008,002A);
- TimezoneOffsetFromUTC (0008,0201);
- LongitudinalTemporalInformationModified (0028,0303);
- DCE private core-lab timing tags (0117,1030)–(0117,1035), method (0117,103A) and comments (0117,103B);
- SER/FTV private sequences (0117,1010) and (0117,10B0).

### Modified dates

TCIA's modified-longitudinal-date process shifts calendar dates while preserving within-patient temporal relationships. The value MODIFIED in (0028,0303) is the explicit indication that dates were shifted. These shifted dates are not real calendar dates. When the flag is present and consistent, use only within-patient intervals (e.g. days from T0 to T1) and keep pseudodates private. If the flag is missing, mixed or inconsistent for a patient/series, the timeline builder fails closed rather than inferring an offset. It preserves date-only precision: a date-only interval is reported at day precision and is not silently upgraded to an exact time-of-day interval.

The DCE private timing tags describe **within-scan** contrast timing (phase starts, injection time, effective delay estimates and SER timing indices). They are not the elapsed time between treatment visits. The DICOM guide cautions that some core-lab injection timing values are best estimates; preserve the AUTO/MANUAL method and comments and do not treat every timing field as exact.

### Trial-stage interpretation

The collection describes the nominal MRI stages as:
- T0: pretreatment (may include a separately consented test/retest scan);
- T1: early treatment, after about three weeks of paclitaxel with or without an experimental agent;
- T2: between paclitaxel and anthracycline/cyclophosphamide (AC);
- T3: after four AC cycles and before surgery.

These are stage definitions, not fixed day intervals. Use actual DICOM dates for patient-specific elapsed intervals. Do not assign 21 days to T0–T1 or assume equal T1–T2/T2–T3 intervals.

## 4. Inventory and timeline commands

Install optional dependencies:

    pip install -e ".[imaging,science,dev]"

Build a local series inventory (private CSV includes IDs, UIDs, shifted dates and filesystem paths):

    python scripts/inventory_ispy2_dicoms.py \
      --dicom-root /private/ispy2_dicom \
      --inventory-csv /private/ispy2_series_inventory.csv \
      --summary-json /private/ispy2_dicom_inventory_summary.json

Build a patient/visit table of valid relative times:

    python scripts/build_ispy2_visit_timeline.py \
      --inventory-csv /private/ispy2_series_inventory.csv \
      --timeline-csv /private/ispy2_visit_timeline.csv \
      --summary-json /private/ispy2_timeline_aggregate.json

The inventory and patient/visit CSVs are patient-level files. Do not commit them. Only aggregate summary JSON is eligible for repository inclusion after a privacy review.

For a human-reviewed series selection, fill templates/ispy2_series_selection_template.csv and validate it against the private inventory. The manifest builder emits a minimal SeriesInstanceUID CSV accepted by TCIA Data Retriever:

    python scripts/prepare_ispy2_series_manifest.py \
      --inventory-csv /private/ispy2_series_inventory.csv \
      --selection-csv /private/ispy2_series_selection.csv \
      --manifest-csv /private/ispy2_selected_series.csv \
      --summary-json /private/ispy2_selected_series_summary.json

Only set include=1 after visually checking the patient, visit and series description. This tool does not infer which image is the correct anatomical/DCE series.

## 5. Anatomical registration and DCE-map resampling

Select a single appropriate 3D anatomical series for each visit (usually a consistent T2-weighted acquisition after checking sequence/protocol metadata). T2-weighted reconstructions can have multiple series; do not choose based on filename order alone. Do not pass a multi-phase DCE acquisition to the rigid-registration script.

Register a follow-up anatomical volume to the T0 reference:

    python scripts/register_ispy2_dicom_series.py \
      --fixed-dicom-dir /private/T0_series \
      --fixed-series-uid T0_SERIES_UID \
      --moving-dicom-dir /private/T1_series \
      --moving-series-uid T1_SERIES_UID \
      --output-image /private/registered_T1_to_T0.nii.gz \
      --output-transform /private/T1_to_T0.tfm \
      --qc-json /private/T1_registration_qc.json

The current implementation is a conservative 3D rigid Euler registration using Mattes mutual information and a multi-resolution schedule. It is a first-pass alignment, not a validated breast-motion solution. Breast positioning, compression, deformation and post-treatment shape changes may not be captured by a rigid transform. Visual overlays, anatomical landmarks and tumour-mask Dice (when independent masks exist) must be reviewed before use. Do not use numerical registration metrics alone as proof of anatomical correctness.

The same-visit DCE PE/SER maps and analysis mask may be resampled with the anatomical transform only when each derived map's DICOM FrameOfReferenceUID has been verified to match that visit's anatomical reference. The resampling helper requires those UID pairs explicitly:

    python scripts/resample_ispy2_dce_map.py \
      --fixed-map-reference /private/T0_SER.nii.gz \
      --moving-map /private/T1_SER.nii.gz \
      --transform /private/T1_to_T0.tfm \
      --output /private/T1_SER_registered.nii.gz \
      --qc-json /private/T1_SER_registration_qc.json \
      --fixed-map-frame-uid T0_SER_FRAME_UID \
      --fixed-anatomical-frame-uid T0_T2_FRAME_UID \
      --moving-map-frame-uid T1_SER_FRAME_UID \
      --moving-anatomical-frame-uid T1_T2_FRAME_UID \
      --kind map

Use --kind label for already decoded categorical masks so nearest-neighbour interpolation is used. The equality check is necessary but not sufficient: confirm physical geometry and inspect overlays. If the DCE map and same-visit anatomy do not share a FrameOfReferenceUID, estimate and validate an additional same-visit map-to-anatomy transform before applying the inter-visit transform. Do not warp masks or maps merely because array dimensions match.

The initial registration is rigid only. A later deformable stage may be considered after quantitative validation on anatomically meaningful structures, but its regularisation and Jacobian/volume-distortion behaviour must be assessed; a deformable field must not silently change the FTV measurement.

## 6. Measurement-side DCE-MRI observation operator

The TCIA DCE derived package includes:
- PEearly: percent enhancement at an effective post-contrast time around 120–150 seconds;
- PElate: percent enhancement around 450 seconds;
- SER: PEearly / PElate, with DICOM rescale slope/intercept required to recover the original SER range;
- a bit-encoded FTV analysis mask;
- private SER/FTV analysis parameters and FTV result sequences.

The implemented measurement-side operator is, given correctly aligned PE/SER maps, the official eligibility mask, voxel spacing and extracted patient-specific thresholds:

\[
\widehat{FTV} =
\Delta v_{\mathrm{cc}}\sum_{x\in\Omega}
\mathbf 1\{x\in M_{\mathrm{pre}}\}
\mathbf 1\{x\in M_{\mathrm{PE}}\}
\mathbf 1\{x\in M_{\mathrm{VOI}}\}
\mathbf 1\{x\notin M_{\mathrm{OMIT}}\}
\mathbf 1\{PE_{\mathrm{early}}(x)\ge PE_{\mathrm{threshold}}\}
\mathbf 1\{SER_{\min}<SER(x)\le SER_{\max}\}.
\]

The maximum SER threshold is included only when defined by the source analysis parameters. SER minimum for the PE-only FTV and the SER-restricted FTV can differ. The PE threshold is nominally 70% in the collection guide, but may be patient-specific and reused at subsequent visits; extract the actual threshold rather than hard-coding 70%.

The DICOM SER map must be rescaled using (0028,1053) Rescale Slope and (0028,1052) Rescale Intercept before thresholding. The DICOM private FTV result sequence (0117,10B0) stores SER minimum, SER maximum, voxel count, volume (cc) and label in child tags (0117,10B1)–(0117,10B5). Recalculate FTV from maps/masks and compare it to the private result sequence as a measurement-QC check.

**Mask caution:** the TCIA analysis mask is bit-encoded and encodes multiple processing layers (precontrast background, PE threshold, VOI and OMIT). It is not an ordinary binary mask. The code's decode_ftv_bitmask function requires explicit, distinct bit assignments supplied by the researcher from the official [Analysis mask files description](https://wiki.cancerimagingarchive.net/download/attachments/50135447/Analysis%20mask%20files%20description.v20211020.docx?api=v2); the repository intentionally does not guess those values. Do not cast the raw mask using mask > 0. The mapping still requires retrieval and verification against that source document before applying the decoder to real cases.

The DICOM FTV computed by UCSF's derived-object pipeline is not guaranteed to be identical to the trial's primary Hologic/AEGIS FTV. Treat it as a related measurement and compare the derived results with the workbook's volume fields before asserting unit or method equivalence. The workbook's VOLUME_TUM_BLU_V1x columns remain unit-unverified.

## 7. Latent PDE-to-MRI observation bridge

The current code implements the measurement-side map/mask-to-FTV operator H_DCE. It does **not** yet implement a calibrated map from latent ONCO-DYNAMICS-X fields to PE/SER. The missing forward model should be staged:

\[
u(x,t)=\big[C,I,F,O,N,D,E,V\big](x,t)
\quad\longrightarrow\quad
\eta(x,t)=g_\theta\big(u(x,t),w(x,t)\big)
\quad\longrightarrow\quad
\widehat{S}(x,\tau)
\quad\longrightarrow\quad
(\widehat{PE}_{\mathrm{early}},\widehat{PE}_{\mathrm{late}},\widehat{SER})
\quad\longrightarrow\quad
\widehat{FTV}.
\]

Here w denotes acquisition- and contrast-kinetic variables; eta must represent or parameterise tissue delivery, vascular volume/permeability and extracellular contrast kinetics (for example Ktrans, ve, vp, if a suitable DCE pharmacokinetic model and arterial input function can be supported). These quantities cannot be read directly from tumour cell density. The MRI signal step must use the actual scanner sequence, baseline T1/flip/TR/TE/rescale metadata and within-scan timing when available. The current PE/SER maps can serve as measured targets; the corresponding g_theta and contrast-signal mapping remain to be selected, calibrated and independently evaluated.

A reasonable next calibration objective is a patient-held-out, voxel-aware loss comparing predicted versus observed PEearly, PElate and SER, plus a secondary FTV agreement term against the DICOM FTV result sequence. Use only training patients to select thresholds, smooth approximations, observation-model coefficients and hyperparameters. All scans from one participant must remain in one partition. Include uncertainty from registration, the bitmask/VOI, signal rescaling and FTV segmentation. Do not describe MRI enhancement volume as a direct viable-cell count.

## 8. Treatment-exposure alignment

The 985-row clinical workbook contains an arm label but not actual administration dates/doses. To estimate exposure-at-scan, obtain a verified event-level administration table from an authorised trial source, with actual administered (not planned) events, drug, dose, units, route, date/time and source provenance.

TCIA shifted DICOM dates cannot be joined to original clinical calendar dates without a documented patient-specific date-shift mapping. The current alignment script therefore refuses ordinary calendar dates and mixed date flags. Populate templates/ispy2_treatment_administration_template.csv only with events already transformed onto the same TCIA modified longitudinal date axis by an authorised, documented process. Days-since-dose and cumulative administered dose are descriptive exposure variables, not drug concentration, AUC, PK/PD or drug-effect coefficients. If required date mapping or actual administration data are unavailable, leave treatment exposure unaligned and do not infer it from Arm.

## 9. Patient-level validation and modelling entry gate

Before fitting the coupled reaction–diffusion model:
1. Verify DICOM-to-workbook subject matching and coverage, with exact one-to-one/ambiguous-match counts.
2. Verify T0–T3 study selection, including optional test/retest scans and multiple studies at one visit.
3. Validate scan timestamp precision and longitudinal-date flags; calculate only supported within-person intervals.
4. Confirm same-visit and across-visit spatial registration using visual and quantitative QC.
5. Verify DCE map slope/intercept, PE timing, SER timing, patient-specific PE/SER thresholds, bitmask layers, voxel geometry and FTV private results.
6. Obtain and align actual treatment administration events, or predefine a model that does not use drug-level exposure.
7. Define the PDE-to-DCE forward observation model; test parameter identifiability on synthetic data and then on training patients.
8. Lock patient-level train/validation/test partitions before fitting; all visits and spatial regions from one participant stay in one partition.
9. Report calibration, prediction error, uncertainty, missingness and excluded cases. Do not report diffusion, proliferation, immune-killing or drug-effect coefficients unless the observation system and identifiability support those interpretations.

## 10. Reproduction and privacy

The DICOM inventory, patient timeline, selected series table, registration outputs, transformations, maps, masks, treatment event tables and participant-level QC should remain in private, approved local/institutional storage. Do not commit patient IDs, raw paths, DICOM UIDs, shifted dates, pixel arrays, registration transforms or the private patient-level workbook.

Public aggregate summaries may be committed only after a privacy review. Current repository tests exercise date parsing, ID/visit parsing, mask/threshold logic, patient-level guards and selected-series manifest checks; real-world registration and DCE output have not been validated because no DICOM images were available locally.

## Source documentation

- [TCIA I-SPY 2 collection and DCE derived objects](https://www.cancerimagingarchive.net/collection/ispy2/)
- [TCIA ACRIN 6698/I-SPY2 imaging protocol, series and mask documentation](https://wiki.cancerimagingarchive.net/pages/viewpage.action?pageId=50135447)
- [ACRIN 6698/I-SPY2 DICOM data descriptions and private tags](https://www.cancerimagingarchive.net/wp-content/uploads/ACRIN-6698-ISPY2-DWI-and-DCE-MRI-Data-Descriptions_20210520.pdf)
- [TCIA longitudinal date de-identification explanation](https://www.cancerimagingarchive.net/submission-and-de-identification-overview/)
- [TCIA Data Retriever source and current CLI usage](https://github.com/TCIA/data-retriever)
- [TCIA API v4 guide](https://wiki.cancerimagingarchive.net/pages/viewpage.action?pageId=70230731)
