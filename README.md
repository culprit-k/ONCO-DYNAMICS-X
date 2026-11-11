# ONCO-DYNAMICS-X

**Multiscale mathematical and AI-assisted breast-cancer research prototype**

ONCO-DYNAMICS-X is a research and university-project prototype that combines dimensionless reaction–diffusion equations, a stochastic cell-agent model, spatial graph features, optional machine-learning architectures, cohort analysis, and an experimental hardware console.

> **Research boundary:** this repository is not a medical device or a clinically validated diagnostic, prognostic, screening, treatment-response, or treatment-selection system. Numerical parameters are provisional and dimensionless unless explicitly documented otherwise. Synthetic or retrospective cohort results are not patient-specific evidence.

## Current project state

| Area | Implemented | Not yet established |
|---|---|---|
| 2D/3D mathematical core | Eight-field reaction–diffusion solver, finite-difference Laplacian, no-flux boundaries, explicit stability guard | Biological parameter calibration, identifiability and clinical validity |
| Cell-agent model | Stochastic cancer, immune and fibroblast agents in 2D/3D | Validation against registered tissue imaging |
| PDE–agent coupling | Operator-split two-way field/agent exchange with explicit occupancy and coupling controls | Mass/number consistency, deposition-kernel convergence and biological calibration |
| AI architectures | PINN, FNO and spatial GNN trained and tested on synthetic ground-truth tasks; expression transformer module | Training/evaluation on biological data and independent validation |
| Cohort pipeline | Public GDC/cBioPortal workflows; 16-marker availability audit and exploratory 50-gene PAM50-panel PC model externally evaluated | Probe-level rederivation, deployable single-sample normalization and independent-cohort molecular validation |
| Survival baseline | Age-only and age + PAM50-panel-PC Cox models evaluated in METABRIC with IPCW time-dependent metrics | Clinical utility, prospective validation, stable molecular improvement and covariate-dependent censoring analysis |
| Biological calibration | Scalar exponential/logistic/Gompertz untreated-volume estimator; local CSV runner; aggregate clinical-target audit | Real longitudinal series, spatial observation operator, biological PDE parameter estimates and external temporal validation |
| Dashboard | Streamlit 2D/3D PDE console and optional coupled agent mode | Full experiment-management layer and validated patient-data ingestion |
| Hardware concept | ESP32 firmware, serial protocol, OpenSCAD phantom and build guide | Physical fabrication, wiring verification and bench testing |

## Verified evaluation results

### Age-only external survival baseline

A separate age-only Cox model was fitted in the TCGA-BRCA PanCancer Atlas cohort and evaluated in METABRIC. Stage was deliberately excluded because TCGA's patient-level AJCC pathologic stage and METABRIC's sample-level generic tumour-stage field have not been shown to be equivalent.

- Development cohort: 1,069 patients, 151 death events.
- External cohort: 1,979 patients, 1,143 death events.
- METABRIC Harrell C-index: **0.588** (patient-bootstrap 95% interval 0.568–0.606).
- IPCW time-dependent AUC at 12 / 36 / 60 months: **0.578 / 0.515 / 0.539**.
- IPCW Brier scores at 12 / 36 / 60 months: **0.0150 / 0.1014 / 0.1713**.
- Sparse-grid integrated Brier score: **0.07935**.

This is weak discrimination and should be treated as a reference comparator, not a useful standalone prognostic model. At five years, the average survival prediction overestimated the IPCW observed survival by about 2.2 percentage points. The IPCW analysis uses a marginal censoring model and assumes independent censoring within METABRIC.

See [full survival results](docs/age_only_external_results.md), [metric methodology](docs/survival_validation_metrics.md), and [aggregate metrics JSON](results/survival_validation/age_only_external_metrics_summary.json).

### PAM50-panel molecular external evaluation

A fixed 50-gene panel was evaluated in TCGA-BRCA (1,069 survival-analysis patients) and METABRIC (1,979 patients). All 50 genes had finite values in both gene-level expression profiles, and the workflow completed an outcome-blind feature pipeline: within-cohort per-gene rank normalisation, TCGA-fitted PCA, and an age + five-PC Cox model evaluated in METABRIC.

- Age-only METABRIC Harrell C-index: **0.5882**.
- Age + PAM50-panel-PC Harrell C-index: **0.5992**.
- Paired difference: **+0.0111**, with 95% bootstrap interval **−0.0019 to +0.0241** (includes zero).
- Age + PC time-dependent AUC at 12 / 36 / 60 months: **0.623 / 0.544 / 0.580**.
- At 60 months, average predicted survival exceeded the IPCW observed estimate by about **3.0 percentage points**.

This is an exploratory, transductive analysis, not a deployable PAM50/Prosigna test. METABRIC's full unlabelled expression distribution is used for rank normalisation; sample-type metadata are unknown for the profile-list entries; probe-level array mapping was not rederived; and the paired C-index difference does not establish stable improvement. See [PAM50 external-validation protocol and results](docs/pam50_molecular_validation.md) and [aggregate report](results/molecular_audit/pam50_pc_external_validation.json).

### PDE–agent coupling sensitivity

The bidirectional coupling prototype completed a 40-scenario synthetic 2D/3D sweep (five coupling strengths × four occupancy scales × two grid dimensions). At full coupling, 80–96% of nodes in the three cell-density fields were at the lower bound across the tested occupancy settings, while upper-cap occupancy remained much smaller. This reveals sparse-grid behaviour and strong dependence on arbitrary bridge parameters; it does not identify a biologically optimal setting. See [sweep results](docs/pde_abm_coupling.md#completed-synthetic-coupling-parameter-sweep--11-november-2026) and [aggregate summary](results/coupling_sensitivity/synthetic_coupling_sweep_summary.json).

### Synthetic PDE numerical validation

A GitHub Actions workflow compared explicit Euler against an adaptive SciPy DOP853 method-of-lines reference on small synthetic grids:

- 2D maximum per-field relative L2 difference: **2.10 × 10⁻⁶**.
- 3D maximum per-field relative L2 difference: **2.35 × 10⁻⁷**.
- Single nested-grid comparison maximum difference: about **0.17% (2D)** and **0.12% (3D)**.

These results check numerical consistency for selected synthetic conditions. They do not validate the biological assumptions, parameters, spatial stencil, patient outcomes or clinical utility. Only one coarse/fine pair per dimension was tested; a three-level convergence study remains necessary.

See [numerical validation results](docs/synthetic_pde_validation_results.md) and [aggregate results JSON](results/numerical_validation/synthetic_pde_validation_summary.json).

### Synthetic AI architecture training (not biological validation)

The optional PINN, FNO and spatial GNN implementations have now completed CPU training on separate analytic/synthetic tasks, with held-out evaluation:

- Boundary-aware PINN (1,000 epochs): held-out relative L2 error 0.02774; PDE residual RMS 0.02297; Neumann boundary residual RMS 0.01188; initial-condition relative error 0.03550.
- FNO: held-out relative L2 error 0.01666, compared with 8.1083 for the no-diffusion persistence baseline on the fixed periodic heat task.
- Spatial GNN: relative L2 error 0.06995 on a synthetic neighbourhood-fraction task; test MSE 0.000687 compared with 0.040535 for the constant-mean baseline.

The 11 November workflow passed tests, Ruff lint and all three training runs. These scores are not breast-cancer accuracy and do not validate the biological model. See [synthetic AI benchmark methodology and results](docs/synthetic_ai_benchmarks.md) and [aggregate metrics](results/ai_benchmarks/synthetic_pinn_fno_gnn_training_summary.json).

### Molecular marker availability audit

A fixed 16-marker panel was queried in the public TCGA-BRCA PanCancer Atlas RSEM profile and METABRIC Illumina HT-12 v3 profile. All **16/16 markers** returned finite expression values for all **1,082 TCGA samples** and all **1,980 METABRIC samples**; no missing, non-finite or duplicate sample records were observed in this run. This is data-availability evidence only. TCGA and METABRIC values are on different assay scales and have not been harmonized or compared directly. See [audit methodology/results](docs/molecular_marker_audit.md) and [aggregate marker report](results/molecular_audit/marker_panel_overlap_summary.json).

### Phase 3 — 50-gene molecular external evaluation

A fixed 50-gene PAM50 research-panel feature inventory was evaluated using TCGA-BRCA for development and METABRIC for external testing. The model used cohort-wise, outcome-blind rank normalisation, five principal components fitted in TCGA only, and a penalised Cox model compared with age-only:

- METABRIC age-only C-index: **0.5882** (95% CI 0.5683–0.6065).
- METABRIC age + molecular PCs C-index: **0.5992** (95% CI 0.5798–0.6201).
- Paired C-index difference: **+0.0111** (95% CI −0.0019 to +0.0241); the interval includes zero, so the gain is not conclusive.
- All 50 genes had finite values in both expression profiles: 54,100 TCGA and 99,000 METABRIC sample–gene measurements.

**Sample selection verified:** the final cBioPortal DETAILED projection confirmed that all 1,082 TCGA and 1,980 METABRIC expression samples were `Primary Solid Tumor`, with zero missing patient links and zero duplicate primary specimens. The expression transformation remains transductive and uses each cohort's full unlabelled expression distribution. This is not a clinical PAM50/Prosigna assay, an intrinsic-subtype classifier, or an individual-patient deployment pipeline.

See [methodology](docs/pam50_molecular_validation.md), [external results and limitations](docs/pam50_molecular_validation_results.md), and [aggregate metrics](results/molecular_validation/pam50_tcga_metabric_summary.json).

A separate five-class subtype-transfer benchmark was trained in TCGA and evaluated in METABRIC using the five shared annotation classes. It achieved **81.0% accuracy**, **74.7% balanced accuracy**, and **77.8% macro-F1**, versus 39.9% accuracy for the TCGA majority-class baseline. METABRIC's claudin-low (217) and NC (6) labels were excluded from this five-class task. Normal-like was the weakest class (F1 0.60), and LumB-to-LumA confusion was common. This is annotation transfer using the same historical subtype gene panel, not independent subtype discovery or an official PAM50/Prosigna classifier. See [methodology](docs/subtype_transfer_benchmark.md) and [full subtype results/confusion matrix](docs/subtype_transfer_results.md).

## Architecture

The current continuum state is:

```
u(x,t) = [C, I, F, O, N, D, E, V]
```

where cancer, immune and fibroblast density proxies interact with oxygen, nutrient, drug, immune-evasion and vasculature proxy fields through a dimensionless reaction–diffusion system.

The optional agent-based system uses cancer, immune and fibroblast agents. In the coupled mode:
1. the PDE advances the continuum fields;
2. agents move, divide or die under local field-dependent rules;
3. agent phenotypes are mapped to grid density proxies;
4. these proxies are blended back into the three cell-density fields.

The agent-to-grid occupancy scale and coupling strength are explicit experimental controls. They are not biologically calibrated quantities. The present blend is an operator-splitting prototype, not a proven mass-conserving multiscale formulation.

## Quick start

Install the numerical core:

```bash
pip install -e .
```

Install scientific and development dependencies:

```bash
pip install -e ".[dev,science]"
pytest
ruff check src tests
```

Install dashboard and survival extras as needed:

```bash
pip install -e ".[dashboard,survival]"
streamlit run dashboard/app.py
```

Run a simulation from the command line:

```bash
onco-dynamics --help
```

Run the synthetic numerical validation:

```bash
python scripts/run_numerical_validation.py
```

Run the age-only survival evaluation (requires access to the public cBioPortal endpoints):

```bash
pip install -e ".[dev,survival]"
python scripts/run_survival_validation_metrics.py
```

## Research documentation

- [Phase 2 validation plan and decision gates](docs/phase_2_validation_plan.md)
- [Implementation status](docs/implementation_status.md)
- [Development roadmap](docs/roadmap.md)
- [Mathematical model and assumptions](docs/mathematical_model.md)
- [PDE–agent coupling method and limitations](docs/pde_abm_coupling.md)
- [Synthetic PDE numerical validation results](docs/synthetic_pde_validation_results.md)
- [Age-only external survival results](docs/age_only_external_results.md)
- [Censoring-aware survival metric methodology](docs/survival_validation_metrics.md)
- [Synthetic AI benchmark results](docs/synthetic_ai_benchmarks.md)
- [Molecular marker audit](docs/molecular_marker_audit.md)
- [PAM50-panel molecular external validation](docs/pam50_molecular_validation.md)
- [PAM50-panel molecular external results](docs/pam50_molecular_validation_results.md)
- [Data sources and governance](docs/data_sources.md)
- [Safety and limitations](docs/safety_and_limitations.md)
- [Hardware build guide](docs/hardware_build_guide.md)

## Automated workflows

The repository includes GitHub Actions for:
- full unit tests and Ruff lint checks;
- GDC paired TCGA-BRCA open-data audit and technical specimen-label baseline;
- TCGA/METABRIC cohort and clinical endpoint audits;
- age-only TCGA-to-METABRIC survival evaluation;
- synthetic 2D/3D PDE numerical validation;
- synthetic PINN/FNO/GNN training benchmarks;
- aggregate-only TCGA/METABRIC molecular marker availability audit;
- exploratory 50-gene PAM50-panel external survival evaluation;
- 2D/3D PDE-agent coupling sensitivity sweep.

The workflows publish aggregate reports only. Do not commit controlled-access or identifiable participant-level data, patient/sample identifiers, raw expression files or row-level clinical records to this public repository.

## Key unresolved research work

1. Resolve endpoint time-origin definitions using each original cohort data guide; retain exclusions for non-positive times rather than changing them.
2. Review the METABRIC stage data dictionary and field timepoint. Keep the age-only model as the reference until stage equivalence is established.
3. Establish a predeclared molecular question, fixed gene/signature set and defensible probe/gene mapping between TCGA RNA-seq and METABRIC microarray without using external outcomes to tune it.
4. Evaluate any molecular predictor externally using discrimination, calibration, uncertainty and failure analyses.
5. Extend numerical tests to three or more grid resolutions and more than one spatial discretization.
6. Establish number/mass accounting and sensitivity for the PDE–agent coupling before interpreting it biologically.
7. Train optional PINN/FNO/GNN models on controlled synthetic ground truth before testing them against real spatial observations.
8. Obtain registered spatial histology or spatial-omics data before claiming validation of cell geometry, immune infiltration, oxygen gradients or drug penetration.
9. Fabricate and bench-test the physical phantom and ESP32 console before claiming a completed hardware demonstration.

## Biological calibration

The estimator, measurement requirements, MRI data pathway and parameter-identifiability gates are documented in [the biological calibration protocol](docs/biological_calibration_protocol.md). No patient- or animal-derived parameter estimates have been generated by this new workflow yet.

## Data and clinical disclaimer

Only public/open resources are used by the committed cohort workflows. Controlled-access participant data require the relevant approvals and must stay outside this public repository. No model output in this project is intended to inform an individual's screening, diagnosis, prognosis, treatment selection, drug dose or clinical management.
