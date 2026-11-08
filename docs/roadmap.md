# Implementation roadmap

The complete architecture is retained. Phases are dependency-aware implementation order, not feature removal.

| Phase | Workstream | Acceptance evidence | Status |
|---|---|---|---|
| 0 | Repository, docs, package, CI, tests | Reproducible install; passing tests | Initial foundation committed |
| 1 | Coupled PDE/ODE engine | Verified solver, positivity, convergence | Prototype implemented; verification continues |
| 2 | PINNs, Bayesian inference, neural operators | Benchmarked against numerical solver | Planned |
| 3 | Agent-based model and spatial GNN | Reproducible cell-level experiments | Planned |
| 4 | Multi-omics and cohort integration | Cohort audit, sample linkage, provenance | Planned |
| 5 | Optimal control, causal analysis, uncertainty | Declared assumptions and robustness checks | Planned |
| 6 | 3D digital twin and dashboard | Computed fields, interactive display, logs | Initial dashboard scaffold |
| 7 | Breast phantom and hardware console | Safe controls, protocol, integration tests | Design stage |
| 8 | End-to-end research demonstration | Reproducible demo, limitations, validation report | Planned |

## Immediate next tasks
1. Add parameterized initial conditions and an explicit scaling scheme.
2. Add grid/timestep convergence experiments and numerical diagnostics.
3. Compare finite differences with an independent solver.
4. Add spatial-image/mask ingestion with metadata validation.
5. Build a cohort audit pipeline before integrating patient datasets.
6. Add PINNs only after reference-solver behavior is understood.
7. Add Bayesian inference and identifiability analysis.
8. Build ABM and GNN modules with spatial ground truth.
9. Integrate multi-omics with patient/sample provenance.
10. Add real-time 3D visualization, then physical controls.
