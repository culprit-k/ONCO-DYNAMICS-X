# System architecture

## Design principle
The conventional numerical model is the transparent reference implementation. AI components are added as challengers or inference layers and must be benchmarked against it. A more complex model is not automatically more accurate.

## Logical layers
| Layer | Responsibility | Candidate implementation |
|---|---|---|
| Mathematical core | Local reactions, spatial transport, numerical integration | NumPy baseline; later SciPy/FEniCSx/FiPy |
| Physics-informed AI | Approximate PDE solutions and infer latent fields | PyTorch PINNs; neural operators |
| Cellular simulation | Individual cancer, immune and stromal agents | PhysiCell or custom ABM |
| Spatial interaction AI | Cell-neighbourhood message passing | PyTorch Geometric GNN |
| Bayesian calibration | Parameter posteriors and uncertainty | PyMC or NumPyro |
| Multi-omics engine | Molecular features and cohort harmonisation | pandas/scikit-learn, then transformer models |
| Causal analysis | Explicit causal estimands and confounding sensitivity | DoWhy/EconML where assumptions are defensible |
| Optimal control | Optimization of hypothetical inputs | CasADi; simulation-only |
| Visualization | Spatial maps, time traces and 3D scene | Streamlit/Plotly; Three.js or Unity |
| Physical console | Knobs, buttons, status and telemetry | ESP32/Raspberry Pi |
| Data governance | Cohort metadata, permissions and provenance | Versioned manifests and local data storage |

## Data flow
1. A configuration defines domain, initial state, parameters and numerical method.
2. The mathematical engine produces spatial fields and time-series summaries.
3. Calibration estimates parameters only when observations support identifiability.
4. AI models consume compatible observations or simulated fields with full provenance.
5. Visualization renders computed values, uncertainty and limitations.
6. Hardware controls send validated parameter changes to the simulation interface.
7. Experiment metadata, seed, environment and commit SHA are logged.

## Data separation
Each patient, specimen, sample, assay and time point needs a traceable identifier and explicit linkage rule. Modalities from different people must never be treated as paired. Spatial claims require spatially resolved data. Synthetic outputs must be labelled synthetic.

## Boundaries
- Core solver does not depend on dashboard or heavyweight AI libraries.
- AI modules are optional and benchmark against the solver.
- Data engine stores schemas and preprocessing, not restricted raw clinical data.
- Hardware controls simulation inputs only and must not expose a clinical treatment interface.
