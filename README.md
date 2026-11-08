# ONCO-DYNAMICS X

**A multiscale mathematical and AI-assisted digital-twin research platform for breast-cancer biology.**

ONCO-DYNAMICS X links mathematical modelling, computational tumour biology, machine learning, multi-omics analysis, interactive visualization, and a physical demonstration console. Its long-term goal is to explore interactions among tumour cells, immune cells, fibroblasts, oxygen, nutrients, vasculature, and hypothetical therapeutic exposure.

> **Research-use notice:** This educational prototype is not a medical device. It is not validated for screening, diagnosis, prognosis, treatment selection, dose calculation, or patient management. Current simulations use dimensionless illustrative parameters and are not patient-specific predictions or clinical recommendations.

## Current software milestone

The numerical baseline is an eight-field reaction-diffusion model with both 2D slice and 3D voxel-grid solvers, no-flux finite differences, explicit time stepping, a dimension-aware diffusion stability guard, bounded fields, and a command-line demo.

Additional implemented components include:
- Stochastic cell-agent simulation using 2D PDE fields as a local environment.
- Public-only GDC file-metadata query and sample-manifest export, plus cohort-manifest validation and identifier-aware expression preprocessing.
- Local parameter sensitivity and Monte Carlo parameter-range uncertainty.
- Optional PyTorch PINN, spatial GNN, Fourier neural operator and expression-transformer architectures.
- Optional PyMC posterior builder, DoWhy causal wrapper and CasADi mean-field control demonstration.
- Streamlit simulation dashboard with a true 3D voxel volume view, as well as ESP32 serial-control integration.

The values are dimensionless teaching proxies, not calibrated biological concentrations. See [implementation status](docs/implementation_status.md) for the exact boundary between implemented components and remaining work.

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

Run a true 3D voxel simulation from the CLI:

```bash
onco-dynamics simulate --dimensions 3 --grid-size 21 --steps 200 --dt 0.0003
```

The CLI writes `results/simulation_summary.csv`. Launch the optional dashboard with:

```bash
pip install -e ".[dashboard,hardware]"
streamlit run dashboard/app.py
```

## Cohort metadata audit

The optional GDC metadata query uses the official API and restricts results to files marked open-access. It does not download the files.

```bash
onco-dynamics gdc-query \
  --project-id TCGA-BRCA \
  --data-type "Gene Expression Quantification" \
  --workflow-type "STAR - Counts" \
  --max-results 100
```

See [cohort data sources and GDC workflow](docs/data_sources.md). Generated metadata CSVs are local and ignored by Git. Never commit controlled or identifiable participant-level data.

For optional AI and research methods, install the relevant extras in pyproject.toml. For example, `pip install -e ".[pinn,spatial]"` installs PyTorch. Do not install every extra blindly; some frameworks require substantial platform-specific dependencies.

## Mathematical core

```
u(x,y[,z],t) = [C,I,F,O,N,D,E,V]^T
du/dt = div(D grad(u)) + R(u,theta) + B q_D(t)
```

The baseline supports 2D and 3D Cartesian-grid models with constant diagonal diffusion, no-flux boundaries, and explicit Euler updates. Read [the mathematical model](docs/mathematical_model.md) and [assumptions and units](docs/assumptions_and_units.md) before interpreting output.

## Full feature architecture retained

Coupled PDE/ODEs; PINNs; agent-based cancer, immune and stromal cells; spatial GNNs; Bayesian estimation; multi-omics AI and transformers; hypoxia and nutrient dynamics; immune evasion; fibroblast/stroma/ECM; drug diffusion; optimal control; causal analysis; uncertainty and sensitivity; neural operators; interactive 3D digital twin; physical console and breast phantom; multi-cohort integration; external validation; reproducible experiments and pathology correlation.

See [architecture](docs/architecture.md), [roadmap](docs/roadmap.md), [hardware console plan](docs/hardware_console.md), and [deployment instructions](docs/deployment.md). Workstreams are staged for verification, not removed. A planned or optional architecture is not a trained/validated model.

## Potential research datasets

Potential sources include [TCGA-BRCA](https://portal.gdc.cancer.gov/), [METABRIC via cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric), [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/), and [Proteomic Data Commons](https://proteomic.datacommons.cancer.gov/). Each cohort must be audited for endpoints, modality, permissions, identifiers, missingness and compatibility. Do not combine modalities from different patients as if paired. Bulk RNA data alone cannot establish spatial positions, oxygen gradients, or direct cell contact.

## Reproducibility

Record configs, dataset versions, inclusion criteria, preprocessing, random seeds, model checkpoints, environment and commit SHA. Separate synthetic simulation tests from biological validation. Report uncertainty, negative results, limitations and failures.

## License

No open-source license has been selected. Until one is added, ordinary copyright restrictions apply.
