# Phase 2 validation plan and decision log

**Project:** ONCO-DYNAMICS-X  
**Status date:** 11 November 2026  
**Purpose:** move from a broad executable architecture to a defensible, reproducible breast-cancer cohort evaluation.

## Current verified state

The public repository contains the mathematical simulation, optional AI architectures, data-engineering workflows, Streamlit dashboard, tests, and hardware design files. The latest recorded GitHub Actions test run on 8 November 2026 completed successfully.

### Evidence already generated

1. **TCGA-BRCA paired tumour/normal technical benchmark**
   - A recorded pilot artifact contains 20 matched patients / 40 samples.
   - The patient-grouped holdout and five-fold grouped cross-validation both report balanced accuracy and ROC-AUC of 1.00.
   - This is a small specimen-label benchmark. It is not a cancer-screening, diagnostic, prognostic, treatment-response, or tumour-dynamics result.
   - The workflow protocol describes selection of up to 48 matched pairs. Do not conflate that intended/selected maximum with the separate 20-pair pilot metrics artifact. Record the exact run artifact and cohort size for every reported result.

2. **TCGA-BRCA PanCancer Atlas survival cohort audit**
   - 1,084 linked patients/samples; 1,082 RNA-seq expression samples.
   - 1,069 expression-linked patients have positive overall-survival follow-up.
   - The audit found 151 deaths in the complete patient-level endpoint set and documented the 13 non-positive times.

3. **METABRIC survival cohort audit**
   - 2,509 linked patients/samples; 1,980 expression samples.
   - 1,979 expression-linked patients have positive overall-survival follow-up.
   - The audit documented endpoint completeness, non-positive times, and clinical-field coverage.

4. **Exploratory clinical Cox transport**
   - The recorded TCGA-to-METABRIC model used age plus a provisional broad-stage mapping.
   - The reported METABRIC Harrell C-index was 0.616 (bootstrap interval 0.592–0.641); age-only comparator on the stage-complete subset was 0.591 (0.570–0.617).
   - A no-stage-requirement age-only sensitivity analysis retained 1,979 METABRIC patients and reported C-index 0.588 (0.568–0.607).
   - These are exploratory discrimination results, not evidence of clinical utility or a validated prognostic product.

## Scientific decision gates

### Gate A — endpoint definition

**Status: provisional; source-level confirmation required.**

The current pipeline parses cBioPortal OS_STATUS as 1 = deceased/event and 0 = living/censored, and uses OS_MONTHS in months. Before the next locked analysis, archive the exact cohort-specific endpoint definition and time origin from the original TCGA PanCancer Atlas and METABRIC study documentation. Keep the non-positive follow-up records quarantined with explicit counts; never replace zero/negative times with arbitrary positive values.

### Gate B — stage harmonization

**Status: BLOCKED for primary cross-cohort prediction.**

The audit shows that TCGA uses patient-level AJCC pathologic stage values (including IA/IB, IIA/IIB, IIIA/IIIB/IIIC and X), whereas METABRIC exposes a sample-level generic TUMOR_STAGE field with values 0–4. The cBioPortal documentation confirms that custom clinical attributes can be study-specific and may be defined at patient or sample level; similarly named or ordinal-looking fields are not proof of semantic equivalence.

Consequences:
- Do not describe the current stage mapping as validated harmonization.
- Do not interpret the stage-only C-index below 0.5 as evidence that advanced disease is protective.
- Treat age-only as the primary transport comparator until the original METABRIC data dictionary/publication establishes the meaning of category 0, staging system, clinical/pathological status, and assessment timepoint.
- If those definitions cannot be reconciled, do not transport a stage-based model.

### Gate C — molecular cross-platform validation

**Status: BLOCKED pending a predeclared feature and harmonization protocol.**

The TCGA PanCancer Atlas RNA-seq RSEM expression profile and METABRIC Illumina HT-12 v3 microarray log2 intensities are different assays and scales. Do not concatenate the matrices, directly reuse coefficients on raw values, or choose genes/thresholds using METABRIC outcomes.

Before fitting a molecular model:
1. Lock one biological question and one endpoint.
2. Define a fixed, literature-supported gene set or pathway score before external evaluation.
3. Specify stable gene-identifier mapping, duplicate handling, missing-gene rules, and a cross-platform transformation using training data only wherever possible.
4. Run cohort/platform diagnostics and a prespecified sensitivity analysis.
5. Freeze the complete TCGA pipeline before opening METABRIC outcome labels for evaluation.
6. Report patient counts, missingness, calibration, discrimination, uncertainty, and failure cases.

### Gate D — survival model quality

**Status: next implementation milestone.**

The next model evaluation should include more than Harrell's C-index:
- calibration at prespecified clinically interpretable time horizons;
- time-dependent discrimination;
- Brier score / integrated Brier score with an explicitly described censoring method;
- proportional-hazards diagnostics and sensitivity to the assumption;
- bootstrap uncertainty at the patient level;
- comparison with an age-only baseline;
- a complete cohort-flow and exclusion report.

Do not report a calibration or time-dependent metric until it has actually been computed on the locked cohort and its censoring assumptions are documented.

### Gate E — mechanistic digital-twin validation

**Status: not yet validated against patient tissue data.**

The presence of PDE, agent-based, PINN, neural-operator, GNN, Bayesian, causal, and control code does not establish that these components are trained or biologically valid. Bulk expression and clinical survival data cannot by themselves validate spatial cell positions, cell-cell contacts, oxygen gradients, immune infiltration geometry, or drug penetration. Those claims require compatible spatially resolved observations and a defined registration/QC protocol.

## Immediate execution order

1. Preserve the current main history and all raw-data privacy boundaries.
2. Complete source-level OS endpoint review and record the exact references/version used.
3. Resolve or formally abandon stage transport; use age-only as the primary external comparator in the meantime.
4. Implement and unit-test time-dependent survival calibration/discrimination with explicit censoring methodology.
5. Re-run the locked TCGA development / METABRIC external evaluation and publish aggregate-only results.
6. Only after the clinical baseline is stable, design the fixed molecular feature pipeline and independent evaluation.
7. Keep PDE/ABM calibration, spatial validation, and hardware fabrication as separate workstreams with their own acceptance tests.

## Data governance

The repository is public. Never commit controlled-access or identifiable participant-level data, patient/sample identifiers, raw expression files, or row-level clinical records. Publish only aggregate metadata, manifests containing non-identifying public file UUIDs where permitted, code, configuration, and reproducible summary metrics.

## References

- cBioPortal clinical file formats: https://docs.cbioportal.org/file-formats/
- cBioPortal data loading and clinical data conventions: https://docs.cbioportal.org/data-loading/
- TCGA PanCancer Atlas study: https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018
- METABRIC study: https://www.cbioportal.org/study/summary?id=brca_metabric
- NCI Genomic Data Commons: https://portal.gdc.cancer.gov/
