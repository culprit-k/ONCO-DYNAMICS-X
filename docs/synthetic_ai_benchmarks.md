# Synthetic AI architecture benchmarks

**Status:** benchmark code and automated workflow are committed; final training metrics are pending a passing CPU workflow.  
**Purpose:** test whether the optional PINN, Fourier neural operator (FNO), and spatial GNN implementations can be trained against mathematical or synthetic ground truth before any biological data are introduced.

These benchmarks are software/model-development checks. They do not train a breast-cancer predictor and must not be reported as biological or clinical performance.

## PINN: analytic 2D no-flux heat equation

The coordinate-to-field network is trained against an analytic solution of:

```
u_t = D * (u_xx + u_yy)
```

on the square `[-1,1] x [-1,1]`, with zero normal derivative at the boundary. The exact non-negative solution is:

```
h(q,t) = (1 - exp(-D*pi^2*t)*cos(pi*q))/2
u(x,y,t) = h(x,t)*h(y,t)
```

where `D=0.02`. This satisfies the heat equation, is non-negative and has no-flux edges. The model is the repository's `FieldPINN` with a Softplus output. The training objective combines supervised data error and automatic-differentiation PDE residual error:

```
loss = data_MSE + 0.1 * PDE_residual_MSE
```

Training, validation and test coordinate samples use separate fixed random seeds. The report contains held-out relative L2 error, residual RMS, initial untrained-model error and validation checkpoints.

## FNO: exact periodic heat semigroup

The FNO learns the map from an initial random spatial field to its exact periodic diffusion result:

```
u_hat(k,t) = u_hat(k,0) * exp(-D * |k|^2 * T)
```

with `D=0.05`, `T=0.2`, and a 16 x 16 grid. The exact Fourier multiplier generates targets, while the repository's `FourierNeuralOperator2D` is trained with mean squared error. Held-out relative L2 error is compared with the identity/persistence baseline (predicting no diffusion). Training/validation/test fields use independent deterministic seeds.

## Spatial GNN: synthetic neighbourhood statistic

Each graph contains 64 nodes with:
- 2D coordinates uniformly sampled in [-1,1];
- one-hot cancer, immune or fibroblast phenotype labels;
- directed edges for neighbour pairs within Euclidean radius 0.4.

The node target is the fraction of each node's neighbours labelled immune. A two-layer `SpatialCellGNN` is trained to estimate this local statistic on disjoint synthetic graphs. Held-out MSE and relative L2 error are compared with a constant-mean training-set baseline.

This benchmark checks graph message passing and neighbourhood aggregation. The target is an engineered statistic, not a biological immune-recruitment probability, cell-contact measurement, or patient-derived label.

## Reproducibility

The workflow `.github/workflows/ai-benchmarks.yml` installs CPU PyTorch, runs the full test suite and Ruff checks, trains the three models with fixed seeds and publishes aggregate metrics only.

Default run:

```bash
pip install "torch>=2.3" --index-url https://download.pytorch.org/whl/cpu
pip install -e ".[dev]"
python scripts/train_synthetic_ai_benchmarks.py \
  --output results/ai_benchmarks/synthetic_pinn_fno_training.json \
  --pinn-epochs 500 \
  --fno-epochs 120 \
  --gnn-epochs 80 \
  --seed 42
```

No checkpoints or row-level data are committed by the workflow.

## Interpretation and next gate

A successful run demonstrates that each architecture can learn its respective controlled synthetic task. It does not establish:
- cancer-cell segmentation or detection;
- immune infiltration or immune evasion in patient tissue;
- tumour growth prediction or survival prediction;
- drug response or treatment benefit;
- accuracy on GDC, TCGA, METABRIC or spatial-omics data;
- biological validity of the PDE/agent-based model.

The next gate for each architecture is separate:
- **PINN:** compare on more than one PDE, include known boundary/initial conditions, and assess residuals and parameter identifiability.
- **FNO:** vary grid resolution, domain and diffusion coefficients; test generalization outside the training parameter distribution.
- **GNN:** compare against non-graph baselines, vary node count and radius, and only then validate against registered segmentation-derived tissue graphs.

Biological model training requires compatible, permissioned datasets and explicit sample linkage/QC. Bulk expression data alone cannot provide spatial cell coordinates or graph ground truth.


## Completed CPU benchmark — 11 November 2026

**Workflow:** [run 42981975418](https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/42981975418)  
**Evaluated commit:** `743608b1906af75bcd8f6e2619c7ecc2d15720a8`  
**Result artifact:** `synthetic-pinn-fno-training` (aggregate-only; retained by GitHub Actions for 30 days)  
**Committed summary:** `results/ai_benchmarks/synthetic_pinn_fno_gnn_training_summary.json`

| Model | Held-out metric | Result | Comparator |
|---|---|---:|---:|
| PINN | Relative L2 field error | 0.09091 | Untrained initialization: 1.3800 |
| PINN | Held-out PDE residual RMS | 0.05325 | Residual is still non-negligible |
| FNO | Relative L2 operator error | 0.01363 | Persistence/no-diffusion baseline: 8.1083 |
| GNN | Relative L2 node-target error | 0.06995 | Test MSE 0.000687 vs constant-mean MSE 0.04054 |

All runs used CPU and seed 42. PINN trained for 500 epochs, FNO for 120, and GNN for 80.

### Readout

- **FNO:** performed well on this 16×16 periodic heat-semigroup task. Cross-resolution and out-of-distribution generalization remain untested.
- **GNN:** learned the engineered neighbour-fraction statistic on separately generated spatial graphs. It has not learned immune-cell behaviour from tissue data.
- **PINN:** learned the analytic field substantially better than random initialization, but its held-out error and PDE residual remain material. It needs further tuning, including sampling, loss balance and explicit boundary-condition assessment.

The workflow's tests and Ruff lint checks passed, and all three training runs completed. These metrics are synthetic software benchmarks only. They are not cancer accuracy, tumour-growth prediction, immune infiltration, treatment response, prognosis, or biological validation.


## Independent rerun — 11 November 2026

A second CPU run was completed after correcting the Neumann analytic-solution test. [Workflow run 42982733597](https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/42982733597) passed the full test suite, Ruff lint, all three training jobs and artifact publication.

The rerun reproduced the same qualitative results:
- PINN held-out relative L2 0.09091; held-out PDE residual RMS 0.05325.
- FNO held-out relative L2 0.01363.
- Spatial GNN held-out relative L2 0.06995; test MSE 0.000687 versus constant-mean baseline 0.04054.

The aggregate output is also recorded in `results/ai_benchmarks/synthetic_pinn_fno_gnn_summary.json`. The earlier completed run and this rerun are independent executions with the same fixed seed and configuration; small floating-point differences are expected across runner/library execution details. Neither run uses patient or cohort data.


## Boundary-aware PINN retraining — current locked result

**This section supersedes the earlier 500-epoch PINN summary above.** The older values are retained as historical baseline results only.

**Workflow:** [run 42983388997](https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/42983388997)  
**Evaluated commit:** `d2886bae8b3487bc68d728da7099fae32a3a6dbf`  
**Artifact:** `synthetic-pinn-fno-training` (aggregate-only)  
**Current committed summary:** `results/ai_benchmarks/synthetic_pinn_fno_gnn_training_summary.json`

### Revised PINN method

The PINN was increased to 1,000 epochs with a 64-unit, four-hidden-layer network. The loss now includes analytic observations, the PDE residual, the initial condition at t=0, and the Neumann normal derivative at each edge:

```
loss = data_MSE
     + 0.1 * PDE_residual_MSE
     + 0.5 * initial_condition_MSE
     + 0.2 * Neumann_boundary_derivative_MSE
```

The best checkpoint was selected by a validation objective combining field and initial-condition errors. Boundary and PDE residuals were assessed on separate held-out coordinates.

| PINN metric | Boundary-aware result |
|---|---:|
| Training epochs | 1,000 |
| Best epoch | 980 |
| Held-out relative L2 field error | **0.02774** |
| Untrained initialisation relative L2 error | 1.4540 |
| Held-out PDE residual RMS | 0.02297 |
| Held-out Neumann boundary residual RMS | 0.01188 |
| Held-out initial-condition relative L2 error | 0.03550 |
| Best validation field MSE | 0.00009926 |
| Best validation initial-condition MSE | 0.00015141 |

The additional constraints improved field error compared with the earlier 500-epoch baseline (0.09091), but the held-out residuals remain nonzero. This is a better synthetic approximation, not a proof of a generally accurate PINN or a cancer model.

### Updated FNO and GNN values from the same run

| Model | Held-out metric | Result | Comparator |
|---|---|---:|---:|
| FNO | Relative L2 operator error | 0.01666 | Persistence/no-diffusion baseline: 8.1083 |
| GNN | Relative L2 node-target error | 0.06995 | Test MSE 0.000687 vs constant-mean MSE 0.04054 |

The FNO training stream was explicitly seeded independently from PINN training so results do not depend on random-number consumption by a preceding model. Its task remains a fixed 16×16 periodic heat operator. The GNN task remains an engineered geometric neighbour-fraction prediction.

The workflow passed full unit tests, Ruff lint, all three training runs and aggregate artifact publication. These three networks have only been trained on analytic/synthetic tasks. The numbers are not breast-cancer accuracy, spatial immune infiltration, clinical prognosis, treatment response or biological validation.
