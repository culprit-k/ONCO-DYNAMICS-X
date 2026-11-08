# Notebook workflow

Use notebooks for exploration only; move reusable logic into `src/onco_dynamics_x`.

Suggested sequence:
1. `01_cohort_audit.ipynb`: provenance, modalities, identifiers, outcomes and permissions.
2. `02_model_exploration.ipynb`: synthetic dynamics and limiting cases.
3. `03_parameter_calibration.ipynb`: identifiability and parameter uncertainty.
4. `04_validation.ipynb`: held-out and external-cohort evaluation.

Do not commit raw patient-level data or private clinical tables. Clear outputs before committing notebooks that might contain sensitive information.
