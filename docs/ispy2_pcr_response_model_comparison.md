# I-SPY 2 pCR-informed response model comparison

**Status: 12 November 2026.** A first patient-separated modelling comparison has been run on the 384-participant MRI feature cohort linked one-to-one to the separate 985-row clinical workbook. The analysis is aggregate-only. Participant identifiers, row-level observations and the linked workbook remain private.

## Objective and estimand

The immediate objective is to test whether a conventional nonlinear response representation improves (a) prediction of the surgical pCR endpoint using only early/inter-regimen MRI observations and (b) forecasting of the final available MRI response state, relative to a simpler linear response reference.

This is a treatment-stage prediction exercise, not a causal treatment-effect analysis and not an estimate of untreated tumour growth. The primary MRI measure is the workbook volume field \`VOLUME_TUM_BLU_V1x\`. Its physical unit and exact equivalence to the trial's primary FTV implementation are not asserted. It is treated only as a positive functional MRI enhancement-volume proxy, and the use of log ratios makes the response features unit-free if the same measurement scale is used across visits.

## Outcome and measurements

The original 384-participant multi-feature MRI publication codes the binary outcome as **1 = pCR** and **0 = non-pCR**. It defines pCR as absence of residual invasive disease in the breast and axillary lymph nodes after neoadjuvant treatment. The published study is the coding basis used here; retain the source-specific mapping when moving to another cohort or workbook version.

For participant i, define the observed stage-specific log responses:

\[
r_{i1}=\log\!\left(\frac{V_{i,T1}}{V_{i,T0}}\right),\qquad
r_{i2}=\log\!\left(\frac{V_{i,T2}}{V_{i,T0}}\right),\qquad
r_{i3}=\log\!\left(\frac{V_{i,T3}}{V_{i,T0}}\right).
\]

These are contrasts between named visits. The script does **not** assign days, equal-duration intervals, a continuous time coordinate, or per-day kinetic rates. The pCR classifier uses T1 and T2 only; T3 is reserved as the target in a separate MRI forecasting task.

## Conventional response models

Both models use the same locked patient split, the same clinical covariates and the same treatment-arm encoding. The clinical predictors are age, workbook-coded HR and HER2, plus treatment arm as a categorical predictor. The models do not use participant IDs or pCR as an input.

### Linear reference

For pCR, a regularized logistic model uses linear terms in \(r_{i1}\), \(r_{i2}\), age, HR, HER2 and treatment-arm indicators:

\[
\operatorname{logit}\Pr(Y_i=1)=\beta_0+\beta_1r_{i1}+\beta_2r_{i2}+\boldsymbol\gamma^T\mathbf{x}_i.
\]

For the separate T3 forecast, ridge regression predicts \(r_{i3}\) from \(r_{i1},r_{i2}\) and the same baseline clinical/treatment variables.

### Quadratic nonlinear response surface

The nonlinear comparator adds the second-order response terms \(r_{i1}^2\), \(r_{i1}r_{i2}\), and \(r_{i2}^2\). The same quadratic terms are added to the T3 forecast model. This is a low-dimensional phenomenological response surface over the observed visit-stage contrasts, not a mechanistic kinetic law.

All scaling and treatment-arm one-hot encoding are fitted inside training-only pipelines. The pCR classifier never receives T3 values. The T3 forecast never receives pCR. This avoids outcome leakage between the two tasks.

## Patient-level validation

A stratified 80/20 holdout was created by splitting the unique participant identifiers after exact one-to-one linkage. The wide modelling table contains one row per participant, so no participant contributes visits to both partitions.

- Linked participants: 384.
- Training participants: 307 (90 pCR, 217 non-pCR).
- Locked test participants: 77 (23 pCR, 54 non-pCR).
- Patient overlap: 0.
- Split seed: \`20261112\`.

The held-out set is internal to this linked cohort. It is not an external cohort or prospective validation. Uncertainty intervals are paired participant bootstraps of this one locked test set (2,000 resamples); they are conditional on the selected split.

## Locked test-set results

| Task / metric | Linear reference | Quadratic nonlinear |
|---|---:|---:|
| pCR ROC AUC | 0.726 | 0.739 |
| pCR AUC bootstrap 95% CI | 0.598–0.836 | 0.614–0.849 |
| pCR Brier score | 0.191 | 0.189 |
| pCR log loss | 0.562 | 0.559 |
| T3 log-volume-ratio forecast RMSE | 0.997 | 1.001 |
| T3 log-volume-ratio forecast MAE | 0.767 | 0.767 |
| T3 forecast R² | 0.452 | 0.448 |

Paired test-set bootstrap comparisons:

- Quadratic minus linear pCR AUC: **+0.013**; 95% CI **−0.018 to +0.050**.
- Quadratic minus linear T3 forecast RMSE: **+0.004**; 95% CI **−0.014 to +0.020**.

The intervals include zero. This run therefore does **not** establish a reliable advantage for the quadratic response surface. Its pCR discrimination is numerically slightly higher in this split, while its T3 forecast RMSE is slightly worse. Model selection should not be based on these small point differences.

## Coupled reaction–diffusion comparator: current identifiability result

The existing ONCO-DYNAMICS-X continuum model is spatial and coupled: tumour, immune and other state fields interact through reaction terms and spatial transport. The currently linked I-SPY 2 workbooks provide whole-lesion volume and other scalar visit summaries, not registered spatial tumour/microenvironment fields. The workbooks also lack actual acquisition timestamps and regimen exposure timelines, and no validated observation operator currently maps the PDE tumour state to the workbook enhancement-volume measurement.

The diffusion limitation is structural. For a tumour field satisfying

\[
\frac{\partial C}{\partial s}=D\nabla^2 C+R(C,\ldots)
\]

with no-flux boundary conditions, integration over the lesion domain \(\Omega\) gives

\[
\frac{d}{ds}\int_\Omega C\,d\Omega
=D\int_\Omega\nabla^2C\,d\Omega+\int_\Omega R(C,\ldots)\,d\Omega
=D\oint_{\partial\Omega}\nabla C\cdot\mathbf n\,dS+\int_\Omega R(C,\ldots)\,d\Omega
=\int_\Omega R(C,\ldots)\,d\Omega.
\]

Thus the aggregate volume balance removes \(D\). Whole-volume trajectories cannot identify a spatial diffusion coefficient under these conditions. Scalar sphericity, diameter, BPE and pCR summaries do not supply the missing spatial fields or independently identify the PDE transport terms.

The PDE branch is therefore reported as **\`blocked_nonidentifiable\`**, with no fitted diffusion coefficient, no PDE score, and no numerical claim that the conventional model outperforms the coupled model. The required inputs for a defensible spatial calibration are:

1. registered spatial tumour/MRI fields and masks across visits;
2. actual acquisition dates/times;
3. aligned treatment start/stop and exposure timelines;
4. a validated observation operator from latent model state to measured MRI enhancement;
5. spatially resolved microenvironment fields appropriate to the coupled state variables (or a revised, explicitly reduced model).

Even after those inputs are obtained, model identifiability and parameter uncertainty must be evaluated before calling any fitted value a biological diffusion or cell-kinetic coefficient.

## Reproduction

Install the project science and development dependencies, then run this locally against the two private workbooks:

\`\`\`bash
pip install -e ".[science,dev]"
python scripts/develop_ispy2_pcr_response_models.py \\
  --features /private/Multi-feature-MRI-NACT-Data.xlsx \\
  --clinical /private/ISPY2-Imaging-Cohort-1-Clinical-Data.xlsx \\
  --output /private/ispy2_pcr_response_model_comparison.json \\
  --test-size 0.20 \\
  --seed 20261112 \\
  --bootstrap 2000

python -m pytest tests/test_ispy2_pcr_response_models.py
\`\`\`

Only aggregate counts and metrics are written to JSON. The script never writes participant IDs or the linked patient-level table. Do not upload either source workbook or the linked workbook to the public repository.

## Data sources

- [TCIA I-SPY 2 collection](https://www.cancerimagingarchive.net/collection/ispy2/): lists the 384-person multi-feature MRI workbook and the separate 985-case clinical workbook, and describes the four MRI feature families.
- [Li et al. (2020), *Predicting breast cancer response to neoadjuvant treatment using multi-feature MRI: results from the I-SPY 2 TRIAL*](https://doi.org/10.1038/s41523-020-00203-7): source for the 384-person feature study, pCR endpoint/coding and MRI feature definitions.

The dataset collection lists CC BY 4.0 for these workbook resources. Preserve the dataset and publication attribution when redistributing aggregate analyses.
