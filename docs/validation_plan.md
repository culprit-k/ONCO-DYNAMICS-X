# Validation and verification plan

## Software verification
- Unit tests for diffusion, reactions, constraints and state shapes.
- Finite-value and positivity checks; deterministic output for deterministic inputs.
- Limiting cases: zero diffusion, zero source and spatially uniform initial states.
- Timestep refinement: compare dt, dt/2 and dt/4.
- Grid refinement on a common observation scale.
- Independent numerical solver comparison.
- Record config, environment, seed and commit SHA.

## Biological validation
A biological claim requires observations corresponding to the quantity modelled. Oxygen dynamics need appropriate oxygen measurements or defensible proxies; spatial immune interactions need spatially resolved assays or validated histology segmentation. Bulk RNA-seq alone cannot establish spatial geometry.

## Cohort strategy
TCGA-BRCA, METABRIC, GEO and CPTAC/PDC are potential resources, not interchangeable datasets. Audit platforms, sample types, endpoints, treatment history, inclusion criteria and permissions. Harmonize only after documenting what is comparable. Reserve an independent cohort for external evaluation.

## Reporting
Report calibration and uncertainty, not only accuracy. Include confidence/credible intervals, failure cases, subgroup limitations and dataset shift. Prognostic performance is not treatment-predictive performance.
