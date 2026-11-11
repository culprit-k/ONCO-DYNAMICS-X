# Biological calibration and parameter-estimation protocol

**Status (11 November 2026): calibration software and data-feasibility workflow started; no patient- or animal-derived ONCO-DYNAMICS-X PDE coefficient has yet been estimated.**

## 1. What the available cohorts can and cannot identify

The TCGA-BRCA PanCancer Atlas and METABRIC work in this repository is predominantly cross-sectional. The expression cohorts contain 1,082 TCGA RNA-seq samples and 1,980 METABRIC microarray samples, with primary-tumour sample type and patient linkage audited. These data support molecular association and cohort-level phenotype analyses. They do not provide registered serial tumour-cellularity maps, oxygen maps, or repeated tumour-volume observations in the current pipeline.

Consequently, the present molecular and survival results cannot identify:
- a spatial cancer-cell diffusion coefficient;
- oxygen-consumption or oxygen-supply kinetics;
- an immune-cell killing rate;
- an immune-evasion production/decay rate;
- a drug-effect parameter or pharmacodynamic EC50;
- patient-specific tumour growth kinetics.

A survival endpoint is not a repeated tumour-growth measurement. Bulk RNA abundance is not a spatial cell-density or oxygen measurement.

The new workflow `scripts/audit_biological_calibration_targets.py` audits public cBioPortal clinical-attribute metadata for candidate Ki-67, ER, PR, HER2, grade, tumour-size, cellularity and hypoxia fields. It outputs aggregate attribute names and coverage only. A name match is not proof of measurement validity: assay method, unit, specimen identity and temporal ordering must be verified before any fitted model uses that field.

## 2. First estimable quantity: whole-tumour untreated growth kinetics

The first estimator is deliberately narrower than the full PDE. For one tumour with repeated untreated volume measurements, it compares three positive-growth laws:

**Exponential**
`V(t)=V_0 e^{r t}`

**Logistic**
`V(t)=K/[1+(K/V_0-1)e^{-r t}]`

**Gompertz**
`V(t)=K\exp[\log(V_0/K)e^{-a t}]`

Here, time is measured in days, volume in mm³, (V_0) is the observed first volume, (r) or (a) is an effective rate per day, and (K) is a fitted carrying-capacity-like scale. The fit minimizes residuals on log volume, corresponding to a multiplicative measurement-error assumption. AICc compares the candidate curves; the selected model receives a parametric residual bootstrap interval and a local Jacobian/SVD identifiability diagnostic.

**Input requirements:** one subject per series; at least five unique, finite time points; strictly positive volumes; documented untreated interval; time in days; volume in mm³. The current CLI defaults to these CSV columns:

```csv
subject_id,time_days,volume_mm3,treatment_status
subject_A,0,120,untreated
subject_A,3,148,untreated
subject_A,6,181,untreated
subject_A,9,221,untreated
subject_A,12,270,untreated
```

The example rows above are schema examples only, not observed research data. The CLI can filter to a user-specified untreated label. If no treatment column is supplied, the analyst must independently verify that every included observation belongs to an untreated interval.

Run locally with:

```bash
pip install -e ".[science]"
python scripts/fit_longitudinal_growth.py \
  --input /path/to/local_longitudinal_volumes.csv \
  --id-column subject_id \
  --time-column time_days \
  --volume-column volume_mm3 \
  --treatment-column treatment_status \
  --untreated-value untreated \
  --bootstrap 500 \
  --output /path/to/local_growth_calibration.json
```

The output contains pseudonymous series keys, fitted parameter summaries, model-selection statistics, uncertainty intervals and identifiability warnings. Raw observations and source identifiers are not exported. Keep the input and outputs in approved local/institutional storage; do not commit patient-level files to this public repository.

The estimator rejects shrinking series. Those belong in a treatment-response model with a time-aligned exposure function; treating treatment-driven shrinkage as untreated proliferation would bias the growth estimate.

### 2.1 Treatment-exposed net functional-volume kinetics

The separate function **fit_treatment_exposed_net_volume_rate()** fits the model:

    log(V(t)) = alpha + k_net * t + error(t)

For the I-SPY 2 FTV field, the TCIA data dictionary defines the DICOM FTV result in **cc (cm³)**. The local runner therefore defaults to column `ftv_cc` and unit `cc`; confirm the workbook/data-dictionary units before changing this. FTV is enhancement-thresholded functional volume, not the whole anatomical tumour volume or a direct tumour-cell count. [TCIA's FTV data dictionary](https://www.cancerimagingarchive.net/wp-content/uploads/ACRIN-6698-ISPY2-DWI-and-DCE-MRI-Data-Descriptions_20210520.pdf) and [I-SPY 2 collection description](https://www.cancerimagingarchive.net/collection/ispy2/).

The slope k_net is the effective net log-volume change per day. It may be negative (net shrinkage) or positive (net growth). Its 95% interval is based on the log-linear regression standard error and Student-t reference distribution. The cohort summary uses a patient-series bootstrap of the median slope.

Use actual acquisition days, not nominal T0/T1/T2/T3 labels. A regimen/subtype grouping can be supplied for descriptive cohort summaries, but this is not a causal treatment comparison. At least three positive-volume visits are required to estimate within-series uncertainty. Zero or below-detection FTV is censored and is excluded until a proper left-censored observation model is implemented; no arbitrary pseudocount is added.

The expected local CSV schema is:

```csv
subject_id,time_days,ftv_cc,regimen
subject_A,0,12.4,standard
subject_A,21,8.1,standard
subject_A,70,3.9,standard
subject_A,140,1.2,standard
```

The rows are schema examples only, not observed patient data. Run locally:

```bash
pip install -e ".[science]"
python scripts/fit_treatment_net_volume.py \
  --input /path/to/local_ispy2_ftv_long.csv \
  --id-column subject_id \
  --time-column time_days \
  --volume-column ftv_cc \
  --volume-unit cc \
  --group-column regimen \
  --bootstrap 2000 \
  --output /path/to/local_net_volume_rates.json
```

This yields treatment-exposed net rates and group-level descriptive intervals. It does not separate proliferation from therapy-induced death or estimate drug EC50. That requires exact exposure histories and additional biological observations. Keep all source rows and output files in approved local/institutional storage.

## 3. Why a fitted growth rate is not yet `cancer_growth`

The current ONCO-DYNAMICS-X PDE uses dimensionless fields and dimensionless time. Its cancer reaction is:

```
dC/dt = rC*C*(1-C/KC)*gO*gN
        - [kI/(1+rhoE*E)]*I*C
        - etaD*gD*C - deltaC*C
```

The fitted (r) from a whole-tumour volume curve is a **net, whole-tumour effective rate**. It mixes intrinsic proliferation, oxygen/nutrient limitation, cell death, immune effects, treatment and measurement/segmentation effects. The unit conversion (r^*=r_{day}\,T_0) only converts a physical rate to a chosen model-time unit; it does not separate those mechanisms. The function `physical_rate_to_dimensionless_rate` is therefore a unit converter, not a parameter assignment.

A future PDE calibration must specify physical length, time and density scales, define a measurement operator from the simulated fields to observed MRI/volume data, and estimate a small identifiable parameter subset. The other parameters should be fixed from validated external measurements or held under justified priors, with sensitivity to those assumptions reported.

## 4. Evidence needed for each mechanism

| Model parameter family | Required matched observations | Why current cross-sectional bulk expression is insufficient |
|---|---|---|
| Cancer proliferation / net growth | Serial untreated volume or cellularity; Ki-67 or other independent proliferation measure | One expression snapshot does not supply a time derivative or a physical growth unit |
| Cancer diffusion / invasion | Registered serial cellularity maps or spatial tumour-boundary changes with physical voxel spacing | Whole-tumour volume discards spatial gradients |
| Oxygen supply/consumption | Quantitative oxygen-sensitive imaging or validated hypoxia measurements, plus vascular information | Bulk HIF1A/CA9 expression does not measure oxygen concentration or consumption rate |
| Immune killing / evasion | Spatially resolved immune-cell density/state, tumour-cell density and longitudinal outcome under documented treatment/exposure | Bulk CD8A/CD274/CXCL9 expression cannot resolve contact, cell state or killing rate |
| Fibroblast activation/suppression | Matched CAF/fibroblast density and stromal/ECM measurements | A bulk stromal marker score does not estimate a transition rate |
| Drug killing / EC50 / field decay | Time-aligned drug exposure or pharmacokinetics and serial response measurements | Treatment-history flags without timing, regimen and concentration do not identify pharmacodynamics |

## 5. Data pathway selected for biological and spatial calibration

The leading public data source is the **I-SPY 2 Breast Dynamic Contrast Enhanced MRI Trial collection at TCIA**:

- The current TCIA collection page reports 719 participants, DCE-MRI and derived parametric maps, MR/SEG data, follow-up and treatment information, and a multi-feature MRI spreadsheet for a 384-participant subset. The compact multi-feature spreadsheet is approximately 130 KB and is the recommended first-stage dataset. The collection includes up to four MRI time points (pre-treatment T0, early treatment T1, mid-treatment T2 and post-treatment T3), with the early scan around three weeks after starting regimen 1.
- For a complete I-SPY 2 Imaging Cohort 1 replication, TCIA directs researchers to use the combined 985-case manifest (719 I-SPY 2 cases plus 266 ACRIN-6698 cases), rather than treating the 719-case collection alone as that full cohort.
- The TCIA collection lists CC BY 4.0 and requires the collection DOI data citation and compliance with the TCIA Data Usage Policy. Verify the exact version and terms again at download time; keep the DOI and attribution with derived results.
- The compact spreadsheet is the preferred first download. It permits a metadata/schema and feature-coverage audit before any multi-terabyte DICOM retrieval. Do not copy participant-level rows to this public GitHub repository.

The 384-patient **Multi-feature MRI NACT Data** spreadsheet is a compact first-stage source (approximately 130 kB) and includes the published DCE-derived feature set (FTV, sphericity, background parenchymal enhancement and longest diameter). For an initial time-series analysis, use this sheet after verifying its actual columns, units, visit-day fields and missingness. If deriving FTV afresh from the DICOM maps, do not silently mix it with the trial's primary FTV values: TCIA notes that its UCSF-derived DICOM FTV can differ from the I-SPY 2 primary commercial-platform FTV, although the difference is generally small. Choose one measurement source for a given analysis and report it. (See [the published feature paper](https://www.nature.com/articles/s41523-020-00203-7) and [TCIA's collection details](https://www.cancerimagingarchive.net/collection/ispy2/).)

**Important model-design implication:** these I-SPY 2 scans are obtained during neoadjuvant treatment. They are not an untreated growth series. The first analysis should estimate a *treatment-exposed net functional-tumour-volume change rate* by interval, retaining exact visit day, treatment arm/regimen and subtype. It must not call this untreated proliferation, drug EC50, or a PDE cell-killing coefficient. Untreated growth and treatment-induced death cannot be separated from a treated volume trajectory without additional measurements/design assumptions.

A complementary **ACRIN 6698 / I-SPY 2 Breast DWI collection** is relevant for diffusion information: TCIA lists 385 subjects and MR/SEG, clinical, follow-up and molecular-test supporting data. ACRIN 6698 was designed to evaluate quantitative diffusion-weighted MRI for neoadjuvant response; the companion study's diffusion imaging and ADC maps may provide a spatially resolved imaging feature alongside DCE-derived functional tumour volume. ADC is a diffusion metric and only a *proxy* for tissue microstructure/cellularity; it must not be treated as a direct cell-density measurement without histologic validation. The exact eligible cases, visit coverage, map availability, segmentation quality and subject overlap between the DCE and DWI collections must be audited from the current manifests before combining them. The DCE and DWI collections are separate collections, and their reported subject counts must not be added or assumed to be the matched sample count. (TCIA ACRIN 6698 and IDC collection descriptions.) 

A 2025 I-SPY 2 mathematical-model analysis reports three visits (baseline V1, early-treatment V2 and post-course V3) for 91 patients across 10 sites. Its design calibrates on V1/V2 and predicts held-out V3. The paper reports that the imaging data came from TCIA and states that figure-generating data are available on reasonable request; therefore this paper is a strong design/benchmark reference, but its exact 91-person analytic subset and map availability must be verified rather than assumed to be identical to the full public TCIA collection.

Additional method references:
- Jarrett et al., *Nature Protocols* (2021), describes serial quantitative MRI acquisition, registration, segmentation, perfusion/diffusion processing, spatial-model calibration and prediction. At least two registered visits are required; a single preprocessed example patient is a software pipeline check, not cohort validation.
- Wang et al., *Journal of Mathematical Biology* (2024), applies sequential Monte Carlo approximate Bayesian computation to in-vivo tumour-volume time series, including mouse breast-cancer xenografts. This supports a likelihood-free uncertainty framework for stochastic models, but its Voronoi ABM and assumptions are not the ONCO-DYNAMICS-X PDE/ABM and its estimates must not be transferred.

External references:
- TCIA I-SPY 2 DCE collection and current access details: https://www.cancerimagingarchive.net/collection/ispy2/
- TCIA ACRIN 6698 / I-SPY 2 Breast DWI collection: https://www.cancerimagingarchive.net/collection/acrin-6698/ (DOI: 10.7937/tcia.kk02-6d95)
- Imaging Data Commons ACRIN-6698 dataset record: https://portal.imaging.datacommons.cancer.gov/collections/acrin_6698/
- TCIA data usage policy: https://www.cancerimagingarchive.net/data-usage-policies-and-restrictions/
- Jarrett et al. 2021 protocol: https://www.nature.com/articles/s41596-021-00617-y
- I-SPY 2 model generalizability analysis (2025): https://doi.org/10.1158/1078-0432.CCR-25-0668
- Tumour-growth ABM calibration using in-vivo volume data and approximate Bayesian computation (2024): https://doi.org/10.1007/s00285-024-02045-4

The cited MRI and ABM models use different state variables, observation operators and governing equations. Their value here is the measurement design, uncertainty treatment and held-out-prediction strategy—not plug-in coefficients or assumed mathematical equivalence.

## 6. Calibration sequence and acceptance gates

1. **Measurement audit:** source/version, access, patient/specimen/timepoint linkage, assay units, treatment timing, inclusion/exclusion and missingness.
2. **Observation operator:** define exactly how model state maps to observed volume, cellularity, oxygenation or cell counts.
3. **Low-dimensional fitting:** fit one or two parameters at a time; keep physically informed nuisance parameters fixed initially.
4. **Identifiability:** inspect Jacobian singular values/condition number, profile likelihood or posterior correlations; demonstrate recovery on synthetic data with known parameters.
5. **Uncertainty:** use bootstrap or Bayesian posterior intervals; do not report point estimates alone.
6. **Temporal holdout:** calibrate on earlier visits only, predict a later visit without refitting, and report error against measurement repeatability.
7. **Patient-level validation:** report performance across subjects, not just a representative example; stratify failure modes and quantify uncertainty.
8. **PDE integration:** only after unit conversion and the observation operator are validated; preserve an uncalibrated baseline for comparison.

### Parameters that must not be fitted simultaneously at the start

The current reaction model contains many confounded terms. For example, observed cancer-field change alone cannot reliably separate `cancer_growth`, `baseline_cancer_loss`, `immune_killing`, oxygen/nutrient limitation and drug killing. Oxygen dynamics cannot separate consumption from supply using a single oxygen snapshot. Drug killing and EC50 cannot be identified without varying/observed exposure. Start with a low-dimensional subset and add parameters only when independent measurements make them identifiable.

## 7. Current implementation and honest status

Implemented:
- exponential/logistic/Gompertz comparison for untreated longitudinal volume data, with log-scale fitting, AICc/BIC, residual bootstrap intervals and Jacobian/SVD identifiability diagnostics;
- a separate treatment-exposed log-volume slope estimator for positive FTV measurements, with a within-series regression interval and subject-series bootstrap summaries;
- a low-dimensional PDE parameter estimator for one to three positive dimensionless reaction coefficients, using registered field snapshots, fixed remaining parameters, weighted nonlinear least squares and a local Jacobian rank/conditioning gate;
- conditional local-linear covariance for the PDE fit, withheld when the fitted subset is rank deficient or ill-conditioned;
- local CSV runners that retain no raw observations or source identifiers in exported summaries;
- an aggregate-only TCGA/METABRIC biomarker-availability audit;
- synthetic parameter-recovery and non-identifiability tests in the dedicated Python 3.11/3.12 CI workflow.

Not yet established:
- no genuine patient or animal series has been fitted by the new ONCO-DYNAMICS-X estimators;
- no patient-derived biological parameter estimate is claimed;
- no MRI/ADC/segmentation fields have been registered to the ONCO-DYNAMICS-X grid or passed through a validated image-to-model observation operator;
- no spatial PDE diffusion, oxygen, immune, stromal or treatment-effect coefficient is biologically calibrated;
- no MRI maps or Ki-67/IHC measurements have been linked to the molecular features.

The next data-dependent gate is to download the compact I-SPY 2 multi-feature MRI workbook and the relevant DCE/DWI/SEG manifests from TCIA into approved local/institutional storage; audit visit coverage, exact case overlap, acquisition days, treatment regimen, FTV availability and censoring; then apply the net-volume estimator to the treated visits with a later visit held out. Spatial PDE fitting must wait for registered quantitative maps, image-to-model units and an independently justified observation operator. Synthetic recovery validates implementation only, not biological accuracy.

## 8. Public clinical-attribute feasibility audit — 11 November 2026

The aggregate-only cBioPortal audit completed successfully in [workflow run 43000949762](https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/43000949762). Matching is restricted to explicit field-name patterns; the initial overly broad ER-status rule was corrected after it incorrectly matched a general cancer-status attribute. The current report is stored as a 30-day aggregate-only workflow artifact and contains no row-level clinical values or identifiers.

| Candidate measurement | TCGA-BRCA PanCancer Atlas | METABRIC |
|---|---|---|
| Ki-67 / proliferation index | No matching field identified | No matching field identified |
| ER status | No verified ER IHC/status field in this metadata audit | ER_IHC patient-level: 2,426 entities; ER_STATUS sample-level: 2,469 |
| PR status | No matching field identified | PR_STATUS sample-level: 1,980 |
| HER2 status | No matching field identified | HER2_SNP6 patient-level: 1,980; HER2_STATUS sample-level: 1,980 |
| Histologic grade | GRADE sample attribute present but 0 entities with values in coverage audit | GRADE sample-level: 2,388 |
| Tumour size | No matching field identified | TUMOR_SIZE sample-level: 2,360 |
| Cellularity / hypoxia | Buffa, Ragnum and Winter hypoxia scores: 1,079 patient entities each | CELLULARITY patient-level: 1,917 |

These are **metadata coverage counts across each full study**, not confirmed matched counts within the expression sample list, and they are not proof that all measurements are available at the same specimen/timepoint. The ER/PR/HER2 fields need source-dictionary review (assay and coding); grade and tumour-size fields need definition and unit review; the METABRIC cellularity field needs a source-level definition and spatial interpretation check. The TCGA hypoxia scores are derived scores, not direct oxygen measurements. A Ki-67 field was not found by the current attribute-name search.

**Decision:** METABRIC is the more promising candidate for exploratory cross-sectional association between molecular features and pathology annotations (ER/PR/HER2/grade/size/cellularity), after sample-to-patient matching and field-definition review. It still does not provide the repeated spatial measurements required to estimate PDE diffusion, oxygen consumption, immune killing or treatment-response kinetics. The TCGA hypoxia scores may support a separately specified cross-sectional hypothesis test, not an oxygen-transport parameter fit.

The next metadata gate is to calculate exact overlap of each candidate field with the selected expression samples and confirm the field definitions from the original study documentation. The next biological-kinetics gate remains an approved longitudinal MRI/volume dataset with timepoint registration and an early-visit calibration / later-visit prediction split.


## 9. Low-dimensional PDE parameter-estimation engine — 11 November 2026

A separate PDE-fitting interface is now implemented in `src/onco_dynamics_x/calibration/pde_parameter_estimation.py`.

It accepts:
- a registered initial `TumourState`;
- observed time points in dimensionless model-time units, starting at zero;
- one or more observed model fields on the same grid and at those time points;
- explicit observation standard deviations for weighted residuals;
- one to three selected positive reaction parameters with explicit bounds;
- fixed values for every other model parameter, grid spacing and solver timestep.

The fitter uses bounded nonlinear least squares in log-parameter space, a fixed-initial-state observation operator, and the existing PDE solver. It reports estimates, weighted residual error, Jacobian singular values/rank/condition number, and a local linearised covariance only when the fitted subset passes the identifiability gate. If the Jacobian is rank deficient or ill-conditioned, covariance is withheld and the fit is flagged. Bounds and observation-noise scales are analyst inputs, not learned biological facts.

**Critical limitation:** the observations must already be registered to the PDE grid, use the same field definitions, and be expressed on the model's dimensionless time/field scales. This code does not convert MRI ADC into cell density, DCE signal into oxygen concentration, or a physical image voxel into the model grid. Those measurement operators and scale conversions must be independently justified before real data can enter this fitter.

The first synthetic recovery test estimates a single `cancer_growth` coefficient from synthetic cancer/oxygen field snapshots generated by the same PDE solver. This verifies implementation and numerical recovery only. It is not biological calibration, and the synthetic initial condition is not a reconstructed patient tumour.

Run validation with:

```bash
pytest -q tests/test_pde_parameter_estimation.py
```

The dedicated biological-calibration CI now includes this module and its synthetic recovery tests. No fitted result from synthetic data is written into `configs/pde.yaml` or used to overwrite the provisional baseline parameters.

**Parameter-estimation strategy:** first fit one coefficient with all other parameters fixed; then test a two-parameter subset only if independent observed fields and the Jacobian support separation. Do not jointly fit growth, baseline loss, immune killing, oxygen limitation and drug effect from tumour volume alone. Diffusion parameters require spatially registered observations; oxygen supply/consumption require repeated oxygen-sensitive measurements; drug killing/EC50 require time-aligned exposure and response data.


## 10. METABRIC molecular-module/pathology association pilot — 11 November 2026

The aggregate-only run completed in [workflow 43003477612](https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/43003477612); its JSON artifact is available from that workflow. The run retrieved all 1,980 METABRIC expression-list samples, verified all 1,980 as Primary Solid Tumor, found one primary sample per patient, and confirmed 99,000/99,000 finite panel measurements. All four expression modules had 100% predefined gene coverage.

The four module means were evaluated against five existing pathology annotations using five-fold stratified patient-level out-of-fold predictions. The fixed classifier was StandardScaler plus class-balanced L2 logistic regression (C=1); no hyperparameter tuning was performed. Bootstrap intervals use 500 stratified patient resamples.

| Target annotation | Patients | Accuracy | Balanced accuracy (95% bootstrap CI) | Macro-F1 |
|---|---:|---:|---:|---:|
| ER_STATUS | 1,980 | 0.940 | 0.944 (0.933–0.954) | 0.921 |
| PR_STATUS | 1,980 | 0.784 | 0.782 (0.764–0.798) | 0.782 |
| HER2_STATUS | 1,980 | 0.970 | 0.979 (0.972–0.986) | 0.937 |
| GRADE | 1,892 | 0.569 | 0.577 (0.546–0.604) | 0.509 |
| CELLULARITY | 1,917 | 0.440 | 0.458 (0.430–0.482) | 0.391 |

**Interpretation and restrictions:**
- The strong ER/HER2 association is expected in part because the molecular modules contain direct receptor-associated genes (ESR1/PGR and ERBB2/GRB7). This is cross-sectional marker concordance with cohort annotations, not independent mechanistic validation or evidence of causality.
- PR association is moderate. Grade and especially cellularity are weakly recovered by these four bulk-expression modules. Do not map these scores to proliferation, stromal, oxygen or tumour-growth rates.
- The study labels are existing cBioPortal annotations. Original assay definitions, scoring rules, specimen provenance and source data dictionary still require review.
- This is internal cross-validation within one cohort, not independent external validation. The cohort-wise rank transform uses the full unlabeled METABRIC expression distribution and remains transductive.
- The result does not establish that the expression modules measure functional receptor activity, cellularity, oxygen concentration, proliferation rate or any PDE coefficient.

**Decision:** retain the four modules as exploratory phenotype descriptors. ER/HER2 scores may be considered receptor-associated expression summaries only; grade/cellularity scores are not sufficiently predictive to serve as model inputs. No PDE/ABM parameter was estimated or changed by this pilot.

Next, review the METABRIC source dictionary and assay provenance, then seek an independent cohort with matching receptor/pathology assays. For actual mechanistic parameter estimation, continue with longitudinal I-SPY 2/ACRIN imaging: validate the observation operator and physical scale first, then estimate only identifiable low-dimensional PDE parameters against registered serial observations.


## 10. Cross-sectional METABRIC pathology association — 11 November 2026

The fixed molecular modules were evaluated against available METABRIC clinical annotations using five-fold stratified patient-level out-of-fold predictions. The workflow verified all **1,980/1,980 primary solid tumour expression samples**, with **99,000/99,000 finite panel measurements**, no duplicate expression records, and one primary sample per patient. All four module groups had 100% prespecified gene coverage.

| Pathology label | Patients | Accuracy | Balanced accuracy | Macro-F1 | Majority accuracy |
|---|---:|---:|---:|---:|---:|
| ER status | 1,980 | 94.0% | 94.4% | 92.1% | 76.1% |
| PR status | 1,980 | 78.4% | 78.2% | 78.2% | 52.5% |
| HER2 status | 1,980 | 97.0% | 97.9% | 93.7% | 87.5% |
| Histologic grade | 1,892 | 56.9% | 57.7% | 50.9% | 50.3% |
| Cellularity category | 1,917 | 44.0% | 45.8% | 39.1% | 50.3% |

The JSON contains the full class counts, confusion matrices and patient-bootstrap intervals: [aggregate METABRIC pathology module association](../results/biological_calibration/metabric_pathology_module_association_summary.json). The GitHub workflow completed successfully: [run 43003477612](https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/43003477612).

**Interpretation matters more than the highest accuracy.** The ER/HER2/PR targets overlap directly with the feature modules: ESR1/PGR are included in the estrogen/luminal module, and ERBB2/GRB7 are included in the HER2 module. The high ER/HER2 scores therefore represent expected marker-to-annotation concordance within the same cross-sectional cohort, not independent validation of a mechanism. Grade association is modest. The cellularity classification is below the majority-class accuracy and does not support using these four modules as a cellularity surrogate.

The bootstrap intervals resample the fixed out-of-fold predictions; they are conditional on the trained fold models and do not include all training-sample/model-selection uncertainty. The rank transformation remains transductive. This was not an external-cohort validation, a prospective analysis, or a causal test.

**Calibration decision:** this analysis supports using the expression modules as descriptive phenotype annotations only. It does not calibrate tumour growth, cellular proliferation kinetics, spatial diffusion, oxygen supply/consumption, immune killing, drug killing, or any PDE coefficient. No Ki-67 measurement was available in the current audit. The next biological-rate gate remains longitudinal MRI/FTV and registered quantitative spatial maps with a validated image-to-model observation operator.
