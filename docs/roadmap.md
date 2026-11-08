# Implementation roadmap

The complete architecture is retained. Phases are dependency-aware implementation order, not feature removal.

| Phase | Workstream | Acceptance evidence | Status |
|---|---|---|---|
| 0 | Repository, docs, package, CI, tests | Reproducible install; passing tests | Implemented |
| 1 | Coupled PDE/ODE engine | Verified solver, positivity, convergence | Prototype implemented; refinement diagnostics added |
| 2 | PINNs, Bayesian inference, neural operators | Benchmarked against numerical solver | Planned |
| 3 | Agent-based model and spatial GNN | Reproducible cell-level experiments | Initial NumPy ABM implemented; GNN planned |
| 4 | Multi-omics and cohort integration | Cohort audit, sample linkage, provenance | Manifest audit implemented; modality ingestion planned |
| 5 | Optimal control, causal analysis, uncertainty | Declared assumptions and robustness checks | Local sensitivity implemented; advanced methods planned |
| 6 | 3D digital twin and dashboard | Computed fields, interactive display, logs | 2D solver dashboard and 3D surface view implemented |
| 7 | Breast phantom and hardware console | Safe controls, protocol, integration tests | ESP32 sketch and validated serial parser implemented |
| 8 | End-to-end research demonstration | Reproducible demo, limitations, validation report | Planned |

## Immediate next tasks
1. Add automated grid/timestep convergence experiments and an independent numerical solver comparison.
2. Add spatial image/mask ingestion with metadata validation.
3. Add a cohort audit pipeline and download only permitted public data.
4. Implement PINN and neural-operator challengers; benchmark them against the numerical reference.
5. Add Bayesian calibration and identifiability analysis.
6. Expand the ABM and implement spatial graph message passing.
7. Integrate multi-omics only with patient/sample provenance.
8. Add causal and optimal-control modules with explicit assumptions.
9. Add 3D tissue simulation and then connect the physical console to live parameter control.
