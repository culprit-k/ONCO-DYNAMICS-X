# ONCO-DYNAMICS X

**A multiscale mathematical and AI-assisted digital-twin research platform for breast-cancer biology.**

ONCO-DYNAMICS X links mathematical modelling, computational tumour biology, machine learning, multi-omics analysis, interactive visualization, and a physical demonstration console. Its long-term goal is to explore interactions among tumour cells, immune cells, fibroblasts, oxygen, nutrients, vasculature, and hypothetical therapeutic exposure.

> **Research-use notice:** This educational prototype is not a medical device. It is not validated for screening, diagnosis, prognosis, treatment selection, dose calculation, or patient management. Current simulations use dimensionless illustrative parameters and are not patient-specific predictions or clinical recommendations.

## Current software milestone

The first milestone is a transparent runnable numerical baseline: an eight-field 2D reaction-diffusion model (cancer, immune, fibroblast, oxygen, nutrients, drug, immune-evasion activity, and vasculature), no-flux finite differences, explicit time stepping, a diffusion stability guard, bounded fields, tests, a CLI demo, and an optional dashboard. Values are dimensionless teaching proxies, not calibrated biological concentrations.

## Quick start

Python 3.11+ is recommended.

```bash
git clone https://github.com/culprit-k/ONCO-DYNAMICS-X.git
cd ONCO-DYNAMICS-X
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -e ".[dev]"
pytest
onco-dynamics simulate --grid-size 51 --steps 400 --dt 0.0003
```

The CLI writes `results/simulation_summary.csv`. Launch the optional dashboard with `pip install -e ".[dashboard]" && streamlit run dashboard/app.py`.

## Mathematical core

```
u(x,y,t) = [C,I,F,O,N,D,E,V]^T
du/dt = div(D grad(u)) + R(u,theta) + B q_D(t)
```

The current implementation is a 2D Cartesian-grid prototype with constant diagonal diffusion, no-flux boundaries, and explicit Euler updates. Read [the mathematical model](docs/mathematical_model.md) and [assumptions and units](docs/assumptions_and_units.md) before interpreting output.

## Full feature architecture (retained roadmap)

1. Coupled reaction-diffusion PDEs and ODEs.
2. Physics-informed neural networks (PINNs), benchmarked against conventional numerical solvers.
3. Agent-based modelling of cancer, immune, fibroblast, and endothelial cells.
4. Spatial graph neural networks for neighbourhood interactions.
5. Bayesian parameter estimation, identifiability, and posterior uncertainty.
6. Multi-omics AI (genomic, transcriptomic, proteomic, pathology and spatial data), including transformer representations.
7. Hypoxia, oxygen and nutrient dynamics.
8. Tumour-immune interactions and immune evasion.
9. Fibroblasts, stroma and extracellular-matrix interactions.
10. Drug diffusion and hypothetical therapeutic-effect simulation.
11. Optimal control of abstract simulation inputs only.
12. Causal analysis with explicit estimands and assumptions.
13. Uncertainty quantification and sensitivity analysis.
14. Neural operators as an advanced extension.
15. Interactive 3D digital twin visualization.
16. Physical hardware console and breast anatomy phantom.
17. Multi-cohort integration and independent biological evaluation.
18. Reproducible configs, tests, versioning, provenance and documentation.
19. Educational pathology correlation where data and access permit.

See [architecture](docs/architecture.md), [roadmap](docs/roadmap.md), and [hardware console plan](docs/hardware_console.md). Workstreams are staged for verification, not removed; a planned module is not claimed to be implemented.

## Potential research datasets

Potential sources include [TCGA-BRCA](https://portal.gdc.cancer.gov/), [METABRIC via cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric), [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/), and [Proteomic Data Commons](https://proteomic.datacommons.cancer.gov/). Each cohort must be audited for endpoints, modality, permissions, identifiers, missingness and compatibility. Do not combine modalities from different patients as if paired. Bulk RNA data alone cannot establish spatial positions, oxygen gradients, or direct cell contact.

Never commit controlled, identifiable, or restricted participant-level data. Keep raw datasets in approved local/institutional storage; commit permitted code, schemas and aggregate audit outputs only.

## Reproducibility

Record configs, dataset versions, inclusion criteria, preprocessing, random seeds, model checkpoints, environment and commit SHA. Separate synthetic software tests from biological validation. Report uncertainty, negative results, limitations and failures.

## Status

**Current milestone:** repository foundation and initial 2D numerical baseline. Further solver verification, calibration, external validation, AI and hardware integration remain future work.

## License

No open-source license has been selected. Until one is added, ordinary copyright restrictions apply.
