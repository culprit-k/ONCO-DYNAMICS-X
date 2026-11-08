# Implementation roadmap

All original workstreams are retained. Phases describe implementation dependency, not feature removal.

| Phase | Workstream | Acceptance evidence | Status |
|---|---|---|---|
| 0 | Repository, docs, package, CI and tests | Reproducible install; automated tests | Implemented |
| 1 | Coupled PDE/ODE engine | Solver checks, positivity, timestep/grid convergence, independent solver | Dimensionless 2D baseline and timestep diagnostic implemented; independent-solver study pending |
| 2 | PINNs, Bayesian inference, neural operators | Trained models benchmarked against numerical solver; posterior diagnostics | PINN/FNO architectures and PyMC model builder implemented; scientific training/calibration pending |
| 3 | Agent-based model and spatial GNN | Cell-level reproducible experiments and spatial validation | Initial NumPy ABM and GNN architecture implemented; richer validation/training pending |
| 4 | Multi-omics and cohort integration | Cohort audit, sample linkage, preprocessing provenance, external cohort | Manifest audit, expression preprocessing and transformer architecture implemented; cohort ingestion/validation pending |
| 5 | Optimal control, causal analysis, uncertainty | Declared assumptions, identifiable estimands, robust sensitivity | Parameter sensitivity, range-based Monte Carlo, DoWhy wrapper and toy CasADi OCP implemented; real-data analyses pending |
| 6 | 3D digital twin and dashboard | Computed spatial fields, interaction, uncertainty and experiment logs | Streamlit dashboard with 2D maps and 3D surface of a 2D field; true 3D solver pending |
| 7 | Breast phantom and hardware console | Safe electronics, serial protocol and physical integration test | ESP32 sketch, parser and dashboard serial hook implemented; physical fabrication/bench test pending |
| 8 | End-to-end university demonstration | Reproducible integrated demo and validation report | Pending |

## Immediate next tasks

1. Run the expanded GitHub Actions suite on the current branch and resolve any failures.
2. Add grid-refinement and independent-solver comparison reports.
3. Audit a permitted cohort and build a complete, reproducible sample manifest.
4. Train the optional AI modules on synthetic benchmark tasks before introducing biological data.
5. Evaluate the PINN, FNO and GNN against transparent numerical/statistical baselines.
6. Build a 3D tissue model only after the 2D solver and spatial data flow are verified.
7. Bench-test the ESP32 serial interface locally and then build the physical console/phantom.
