# Age-only TCGA-to-METABRIC survival validation

**Implementation status:** code and workflow committed; external result is pending a successful GitHub Actions run.  
**Date:** 11 November 2026

## Question

How well does a simple, frozen clinical baseline based only on age at diagnosis rank and estimate overall survival in the independent METABRIC cohort after fitting in TCGA-BRCA?

This is a clinical-covariate baseline, not the molecular or spatial ONCO-DYNAMICS-X digital twin. It is included as a reference point against which more complex models must eventually be compared.

## Locked design

- Development cohort: cBioPortal brca_tcga_pan_can_atlas_2018, RNA-seq mRNA sample list.
- External cohort: cBioPortal brca_metabric, mRNA sample list.
- Endpoint fields: OS_MONTHS and OS_STATUS.
- Event convention: 1:DECEASED = 1; 0:LIVING = 0.
- Predictor: age at diagnosis only.
- Age standardization: mean and standard deviation fitted in TCGA and reused unchanged in METABRIC.
- Model: Cox proportional hazards with penalizer 0.1.
- Stage: deliberately excluded because the available TCGA patient-level AJCC pathologic stage and METABRIC sample-level generic tumour-stage field have not been established as semantically equivalent.
- Expression: not used.
- External outcomes: not used for fitting or model tuning.
- Evaluation horizons: 12, 36 and 60 months.
- Bootstrap: 250 patient-resampling replicates, seed 42.
- Report: aggregate counts and metrics only. No patient or sample identifiers, row-level clinical records, or expression matrices are written.

## Endpoint provenance and limits

cBioPortal's clinical file-format documentation specifies OS_STATUS prefixes (1 for an event, 0 for no event) and describes OS_MONTHS as months since initial diagnosis. The cBioPortal FAQ explicitly warns that study-specific endpoint definitions may vary, so the original study guides remain controlling.

For TCGA, the Pan-Cancer Clinical Data Resource paper defines overall survival from diagnosis to death from any cause, with censoring from diagnosis to last contact. For METABRIC, public dataset documentation describes OS months as overall-survival months since initial diagnosis. The project should retain those source references with each result version and recheck them before any clinical interpretation.

Non-positive times are excluded by the existing cohort assembler and counted in aggregate quality-control output. They are never replaced with arbitrary positive times.

## Statistical metrics

### Harrell's concordance index

The report includes Harrell's C-index as a global ranking summary. Higher values indicate better concordance of predicted risk with observed event ordering under the estimator's comparable-pair rules. It is not a calibration measure.

### Cumulative/dynamic time-dependent AUC

At a horizon t:
- Cases are patients with an observed event by t: event = 1 and observed time <= t.
- Controls are patients known to be event-free beyond t: observed time > t.
- Patients censored before or at t without a prior event are not directly classifiable at t and receive no case/control label.

Inverse probability of censoring weights are applied:
- Case weight: 1 / G(T_i-), where G is the Kaplan-Meier estimate of the censoring survival function immediately before the case's event time.
- Control weight: 1 / G(t).

The weighted AUC is the weighted probability that a randomly selected case has a greater predicted risk than a randomly selected control, with half credit for ties. This is a cumulative/dynamic definition; it is not incident/dynamic AUC.

### IPCW Brier score

Let S_i(t) be the model-predicted survival probability. The score is calculated as:

```
Brier(t) = 1/n * [
  sum_cases w_i * (0 - S_i(t))^2
  + sum_controls w_i * (1 - S_i(t))^2
]
```

where the case/control weights are the inverse censoring probabilities above. The denominator is the full evaluation-cohort size n, as in the IPCW Brier estimator; early-censored subjects contribute zero to the observed sum because their survival status at t is unknown.

### Calibration summaries

Ten approximately equal-sized groups are formed by predicted survival among patients whose status at the given horizon can be evaluated. Each group reports:
- mean predicted survival, with IPCW weighting;
- IPCW estimate of observed survival;
- number of evaluable patients.

The report also gives the overall IPCW observed-minus-predicted survival gap. These are descriptive calibration summaries, not a fitted recalibration model.

### Integrated Brier score

The workflow reports a sparse-grid trapezoidal approximation using time zero (Brier score anchored at zero) and 12, 36 and 60 months. It is explicitly not a dense-grid integrated Brier score. A future analysis may calculate a denser curve once follow-up support and horizon stability have been assessed.

## Censoring assumptions and uncertainty

The marginal Kaplan-Meier censoring distribution is estimated within the external evaluation cohort for IPCW correction only. The fitted Cox coefficients and individual survival predictions remain frozen from TCGA. This metric procedure assumes independent censoring in the evaluation cohort; it does not model covariate-dependent censoring. If censoring depends on age or other predictors, a conditional censoring model or sensitivity analysis will be needed.

Bootstrap confidence intervals use patient-level resampling of the external cohort and recompute the censoring weights in each replicate. Replicates that cannot estimate a metric because of absent cases/controls or unstable censoring survival are omitted. The report exposes requested and valid replicate counts. An interval is withheld when fewer than half of the requested replicates or fewer than 30 replicates are valid.

A small censoring-survival probability can make inverse weights unstable. The implementation marks a horizon not estimable when the censoring survival probability falls below 0.05 rather than reporting an unstable metric.

## Reproducibility

Run locally or in GitHub Actions:

```bash
pip install -e ".[dev,survival]"
pytest
python scripts/run_survival_validation_metrics.py \
  --output results/survival_validation/age_only_external_metrics.json \
  --horizons-months 12 36 60 \
  --bootstrap 250 \
  --seed 42 \
  --penalizer 0.1
```

The GitHub Actions workflow is .github/workflows/survival-validation.yml. It runs the repository tests and lint checks, retrieves the public cohort data at run time, and publishes only the aggregate report as a short-lived workflow artifact.

## Interpretation boundary

A successful run would establish that the declared workflow executed and provide an exploratory external evaluation of an age-only Cox baseline. It would not establish clinical utility, prospective generalizability, a treatment effect, treatment response, molecular predictive value, or biological validity of the PDE/agent-based/PINN/FNO/GNN components.

The project must not advance to a molecular claim until cross-platform feature mapping and preprocessing are locked without METABRIC outcome leakage. It must not claim that bulk clinical and expression cohorts validate spatial oxygen gradients, immune-cell geometry, cell-cell contacts, or drug penetration.

## References

- cBioPortal clinical file formats: https://docs.cbioportal.org/file-formats/
- cBioPortal FAQ on OS_STATUS / OS_MONTHS: https://docs.cbioportal.org/user-guide/faq/
- TCGA Pan-Cancer Clinical Data Resource: https://pmc.ncbi.nlm.nih.gov/articles/6066282/
- METABRIC OS field descriptions: https://pmacdasci.github.io/IntroR/vignettes/5_Exploring_Metabric.html
- TCGA PanCancer Atlas: https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018
- METABRIC: https://www.cbioportal.org/study/summary?id=brca_metabric
