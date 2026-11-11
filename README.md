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
| AI architectures | Optional PINN, FNO, spatial GNN and expression transformer modules | Training and independent biological evaluation of these models |
| Cohort pipeline | Public GDC and cBioPortal metadata/data workflows; aggregate-only audits | Molecular cross-platform validation and harmonized feature space |
| Survival baseline | TCGA-trained / METABRIC-tested age-only Cox baseline with censoring-aware evaluation | Clinical utility, prospective validation, molecular predictor and covariate-dependent censoring analysis |
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

### Synthetic PDE numerical validation

A GitHub Actions workflow compared explicit Euler against an adaptive SciPy DOP853 method-of-lines reference on small synthetic grids:

- 2D maximum per-field relative L2 difference: **2.10 × 10⁻⁶**.
- 3D maximum per-field relative L2 difference: **2.35 × 10⁻⁷**.
- Single nested-grid comparison maximum difference: about **0.17% (2D)** and **0.12% (3D)**.

These results check numerical consistency for selected synthetic conditions. They do not validate the biological assumptions, parameters, spatial stencil, patient outcomes or clinical utility. Only one coarse/fine pair per dimension was tested; a three-level convergence study remains necessary.

See [numerical validation results](docs/synthetic_pde_validation_results.md) and [aggregate results JSON](results/numerical_validation/synthetic_pde_validation_summary.json).

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
- [Data sources and governance](docs/data_sources.md)
- [Safety and limitations](docs/safety_and_limitations.md)
- [Hardware build guide](docs/hardware_build_guide.md)

## Automated workflows

The repository includes GitHub Actions for:
- full unit tests and Ruff lint checks;
- GDC paired TCGA-BRCA open-data audit and technical specimen-label baseline;
- TCGA/METABRIC cohort and clinical endpoint audits;
- age-only TCGA-to-METABRIC survival evaluation;
- synthetic 2D/3D PDE numerical validation.

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

## Data and clinical disclaimer

Only public/open resources are used by the committed cohort workflows. Controlled-access participant data require the relevant approvals and must stay outside this public repository. No model output in this project is intended to inform an individual's screening, diagnosis, prognosis, treatment selection, drug dose or clinical management.
