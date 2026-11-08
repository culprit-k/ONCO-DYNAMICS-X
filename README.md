# ONCO-DYNAMICS X

**A multiscale mathematical and AI-assisted digital-twin research platform for breast-cancer biology.**

ONCO-DYNAMICS X links mathematical modelling, computational tumour biology, machine learning, multi-omics analysis, interactive visualization, and a physical demonstration console. Its long-term goal is to explore interactions among tumour cells, immune cells, fibroblasts, oxygen, nutrients, vasculature, and hypothetical therapeutic exposure.

> **Research-use notice:** This educational prototype is not a medical device. It is not validated for screening, diagnosis, prognosis, treatment selection, dose calculation, or patient management. Current simulations use dimensionless illustrative parameters and are not patient-specific predictions or clinical recommendations.

## Current software milestone

The numerical baseline is an eight-field reaction-diffusion model with both 2D slice and 3D voxel-grid solvers, no-flux finite differences, explicit time stepping, a dimension-aware diffusion stability guard, bounded fields, and a command-line demo.

Additional implemented components include:
- Stochastic cell-agent simulation using 2D PDE fields as a local environment.
- Public-only GDC file-metadata query, sample-manifest export, manifest audit, and explicit single-file open-data download with checksum verification.
- Identifier-aware expression preprocessing, local parameter sensitivity, and Monte Carlo parameter-range uncertainty.
- Optional PyTorch PINN, spatial GNN, Fourier neural operator and expression-transformer architectures.
- Optional PyMC posterior builder, DoWhy causal wrapper and CasADi mean-field control demonstration.
- Streamlit dashboard with true 3D voxel-volume rendering and ESP32 serial-control integration.

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

Run a true 3D voxel simulation:

```bash
onco-dynamics simulate --dimensions 3 --grid-size 21 --steps 200 --dt 0.0003
```

## TCGA-BRCA public metadata audit

Query open-access GDC metadata (no file downloads):

```bash
onco-dynamics gdc-query \
  --project-id TCGA-BRCA \
  --data-type "Gene Expression Quantification" \
  --workflow-type "STAR - Counts" \
  --max-results 100
```

After reviewing the local metadata CSV, download one explicitly selected open-access file by its GDC UUID:

```bash
onco-dynamics gdc-download-open \
  --metadata-csv data/metadata/gdc_open_files.csv \
  --file-id YOUR-GDC-FILE-UUID \
  --destination data/raw/gdc
```

The downloader refuses non-open files and verifies the byte count and MD5 checksum before finalizing a file. It does not provide controlled-access data or bypass access requirements. Generated CSVs and raw data are ignored by Git. See [cohort data sources](docs/data_sources.md).

## Dashboard

```bash
pip install -e ".[dashboard,hardware]"
streamlit run dashboard/app.py
```

The dashboard can run locally without hardware. For the ESP32 controls, connect the development board over USB to the same computer running Streamlit and use the detected serial device path. Cloud-hosted dashboards cannot directly access a USB port attached to a separate computer.

## Mathematical core

```
u(x,y[,z],t) = [C,I,F,O,N,D,E,V]^T
du/dt = div(D grad(u)) + R(u,theta) + B q_D(t)
```

The baseline supports 2D and 3D Cartesian-grid models with constant diagonal diffusion, no-flux boundaries, and explicit Euler updates. Read [the mathematical model](docs/mathematical_model.md) and [assumptions and units](docs/assumptions_and_units.md) before interpreting output.

## Full feature architecture retained

Coupled PDE/ODEs; PINNs; agent-based cancer, immune and stromal cells; spatial GNNs; Bayesian estimation; multi-omics AI and transformers; hypoxia and nutrient dynamics; immune evasion; fibroblast/stroma/ECM; drug diffusion; optimal control; causal analysis; uncertainty and sensitivity; neural operators; interactive 3D digital twin; physical console and breast phantom; multi-cohort integration; external validation; reproducible experiments and pathology correlation.

See [architecture](docs/architecture.md), [roadmap](docs/roadmap.md), [hardware console plan](docs/hardware_console.md), and [deployment instructions](docs/deployment.md). Workstreams are staged for verification, not removed. A planned or optional architecture is not a trained/validated model.

## Reproducibility

Record configs, dataset versions, inclusion criteria, preprocessing, random seeds, model checkpoints, environment and commit SHA. Separate synthetic simulation tests from biological validation. Report uncertainty, negative results, limitations and failures.

## License

No open-source license has been selected. Until a license is added, ordinary copyright restrictions apply.
