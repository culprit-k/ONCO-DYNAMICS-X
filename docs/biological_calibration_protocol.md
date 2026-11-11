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

- The current TCIA collection page reports 719 participants, DCE-MRI and derived parametric maps, MR/SEG data, follow-up and treatment information, and a multi-feature MRI spreadsheet for a 384-participant subset. The full image collection is approximately 1.75 TB; the compact feature spreadsheet is approximately 130 KB. The collection includes up to four MRI time points (pre-treatment T0, early treatment T1, mid-treatment T2 and post-treatment T3), with the early scan around three weeks after starting regimen 1.
- For a complete I-SPY 2 Imaging Cohort 1 replication, TCIA directs researchers to use the combined 985-case manifest (719 I-SPY 2 cases plus 266 ACRIN-6698 cases), rather than treating the 719-case collection alone as that full cohort.
- The TCIA collection lists CC BY 4.0 and requires the collection DOI data citation and compliance with the TCIA Data Usage Policy. Verify the exact version and terms again at download time; keep the DOI and attribution with derived results.
- The compact spreadsheet is the preferred first download. It permits a metadata/schema and feature-coverage audit before any multi-terabyte DICOM retrieval. Do not copy participant-level rows to this public GitHub repository.

**Important model-design implication:** these I-SPY 2 scans are obtained during neoadjuvant treatment. They are not an untreated growth series. The first analysis should estimate a *treatment-exposed net functional-tumour-volume change rate* by interval, retaining exact visit day, treatment arm/regimen and subtype. It must not call this untreated proliferation, drug EC50, or a PDE cell-killing coefficient. Untreated growth and treatment-induced death cannot be separated from a treated volume trajectory without additional measurements/design assumptions.

A 2025 I-SPY 2 mathematical-model analysis reports three visits (baseline V1, early-treatment V2 and post-course V3) for 91 patients across 10 sites. Its design calibrates on V1/V2 and predicts held-out V3. The paper reports that the imaging data came from TCIA and states that figure-generating data are available on reasonable request; therefore this paper is a strong design/benchmark reference, but its exact 91-person analytic subset and map availability must be verified rather than assumed to be identical to the full public TCIA collection.

Additional method references:
- Jarrett et al., *Nature Protocols* (2021), describes serial quantitative MRI acquisition, registration, segmentation, perfusion/diffusion processing, spatial-model calibration and prediction. At least two registered visits are required; a single preprocessed example patient is a software pipeline check, not cohort validation.
- Wang et al., *Journal of Mathematical Biology* (2024), applies sequential Monte Carlo approximate Bayesian computation to in-vivo tumour-volume time series, including mouse breast-cancer xenografts. This supports a likelihood-free uncertainty framework for stochastic models, but its Voronoi ABM and assumptions are not the ONCO-DYNAMICS-X PDE/ABM and its estimates must not be transferred.

External references:
- TCIA I-SPY 2 collection and current access details: https://www.cancerimagingarchive.net/collection/ispy2/
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
- exponential/logistic/Gompertz comparison for untreated longitudinal volume data;
- log-scale nonlinear least-squares estimation with positive-rate/capacity parameterization;
- AICc/BIC and Jacobian/SVD local identifiability diagnostics;
- parametric residual bootstrap confidence intervals;
- explicit conversion of physical rates to a model-time unit without silently assigning them to PDE coefficients;
- local CSV CLI with pseudonymous subject keys;
- aggregate-only TCGA/METABRIC calibration-target metadata audit;
- synthetic parameter-recovery unit tests and dedicated CI workflows.

Not yet established:
- no genuine patient or animal volume series has been fitted by this new ONCO-DYNAMICS-X code;
- no new biological parameter estimate is claimed;
- no spatial PDE coefficient is calibrated;
- no therapy-response parameter is calibrated;
- no MRI maps or Ki-67/IHC measurements have been joined to molecular features.

The next data-dependent gate is to run the aggregate clinical-attribute audit, review its matched measurement candidates, then obtain an approved longitudinal MRI/cellularity dataset with at least two calibration timepoints and a later holdout. Until then, the new estimator is validated only for mathematical recovery on synthetic data.


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
