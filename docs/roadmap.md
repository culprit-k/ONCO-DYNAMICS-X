# Implementation roadmap

All original workstreams are retained. Phases describe implementation dependencies, not feature removal.

| Phase | Workstream | Acceptance evidence | Status |
|---|---|---|---|
| 0 | Repository, docs, package, CI and tests | Reproducible install; automated tests | Implemented |
| 1 | Coupled PDE/ODE engine | Solver checks, positivity, timestep/grid convergence, independent solver | Dimensionless 2D/3D baseline and timestep diagnostic implemented; independent-solver/grid study pending |
| 2 | PINNs, Bayesian inference, neural operators | Trained models benchmarked against numerical solver; posterior diagnostics | PINN/FNO architectures and PyMC model builder implemented; scientific training/calibration pending |
| 3 | Agent-based model and spatial GNN | Cell-level reproducible experiments and spatial validation | Initial NumPy 2D ABM and GNN architecture implemented; richer validation/training pending |
| 4 | Multi-omics and cohort integration | Cohort audit, sample linkage, preprocessing provenance, external cohort | GDC open-file metadata query, sample manifest, expression preprocessing and transformer architecture implemented; other cohort integration pending |
| 5 | Optimal control, causal analysis, uncertainty | Declared assumptions, identifiable estimands, robust sensitivity | Parameter sensitivity, range-based Monte Carlo, DoWhy wrapper and toy CasADi OCP implemented; real-data analyses pending |
| 6 | 3D digital twin and dashboard | Computed 3D spatial fields, interaction, uncertainty and experiment logs | True 3D reaction-diffusion voxel solver and volume rendering implemented; full twin integration pending |
| 7 | Breast phantom and hardware console | Safe electronics, serial protocol and physical integration test | ESP32 sketch, parser and dashboard serial hook implemented; physical fabrication/bench test pending |
| 8 | End-to-end university demonstration | Reproducible integrated demo and validation report | Pending |

## Immediate next tasks
1. Run GitHub Actions on the current commit and resolve any failures.
2. Add grid-refinement and independent-solver comparison reports.
3. Audit a permitted breast-cancer cohort and review the metadata manifest before downloading data.
4. Train optional AI modules on synthetic benchmark tasks before introducing biological data.
5. Evaluate PINN, FNO and GNN against transparent numerical/statistical baselines.
6. Add spatial histology/mask workflows and modality-specific QC.
7. Bench-test the ESP32 locally and then build the physical console and phantom.
