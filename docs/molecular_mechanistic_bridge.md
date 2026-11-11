# Molecular-to-mechanistic bridge: specification and limits

## Purpose

The fixed-panel molecular analysis now has an explicit, testable interface to the mathematical model. The first implementation computes four descriptive expression modules from the existing cohort-wise rank-normalised panel:

- **Proliferation-associated expression:** ANLN, BIRC5, CCNB1, CDC20, CENPF, CEP55, KIF2C, MKI67, MYBL2, TYMS and UBE2C.
- **Estrogen/luminal-associated expression:** BAG1, BCL2, ESR1, FOXA1, MAPT, MLPH, NAT1, PGR and SLC39A6.
- **HER2-associated expression:** ERBB2 and GRB7.
- **Basal/epithelial-associated expression:** CDH3, EGFR, FOXC1, KRT14, KRT17, KRT5, MMP11 and TMEM45B.

For each sample and module, the score is the arithmetic mean of the available member-gene normal scores. At least 75% of a module's predefined genes must be present. Inputs must be finite and uniquely gene-labelled. No outcomes, subtype labels or fitted coefficients are used.

These are deliberately called *expression modules*, not validated pathway scores. The groupings are transparent heuristics over a fixed research panel; they are not independently validated gene-set signatures, and mean expression does not establish pathway activation.

## Proposed biological interpretation

| Expression module | Candidate model-facing interpretation | Not established by expression alone |
|---|---|---|
| Proliferation-associated | A prior or covariate for cancer-cell proliferation propensity / ABM phenotype state | A numerical cell-division rate, tumour growth rate or doubling time |
| Estrogen/luminal-associated | Luminal/estrogen-associated phenotype state or cell-state composition | Functional ER signalling, endocrine sensitivity or response |
| HER2-associated | HER2-associated phenotype annotation | HER2 amplification, receptor phosphorylation or pathway activity |
| Basal/epithelial-associated | Basal/epithelial phenotype state or cell-state composition | Invasion rate, metastatic potential or treatment resistance |

The machine-readable policy is exposed by `describe_mechanistic_link_policy()` in `src/onco_dynamics_x/data_engine/molecular_mechanistic_bridge.py`.

## Explicitly not connected yet

The module does **not** alter any PDE or ABM parameter. In particular, it does not map expression to:

- diffusion coefficients;
- oxygen consumption, nutrient uptake or hypoxia kinetics;
- immune killing or immune-evasion rates;
- drug sensitivity, pharmacodynamics or treatment response;
- a patient-specific tumour growth rate.

This separation is intentional. Expression values are not direct measurements of these physical or biological rates. The existing PDE fields are dimensionless proxies and their parameters are not yet calibrated to patient-matched measurements.

## Required calibration before coupling

1. Define an observable target for the proposed parameter (for example, matched proliferation assay, Ki-67, longitudinal volume change, or measured oxygenation).
2. Assemble matched molecular and target measurements with documented specimen identity and compatible time points.
3. Prespecify a mapping model and fit it only in a designated training cohort.
4. Quantify identifiability, parameter uncertainty, sensitivity to preprocessing and uncertainty in the measured target.
5. Freeze the mapping and evaluate it once in an independent cohort.
6. Compare against a no-molecular-input baseline and report calibration, uncertainty and failure cases.

For the proposed links, receptor status should be corroborated by IHC/ISH or copy-number data; ER/PR status by pathology; and proliferation by an independent assay. A subtype label is not a substitute for these direct measurements.

## Processing caveat

The current source expects already rank-normalised expression. The existing cohort-wise transform uses the full unlabeled expression distribution of each cohort and is therefore transductive. It is not a prospective single-sample transform. Scores from different cohorts should not be interpreted as absolute comparable biological units. The frozen-reference transformation and array probe-level METABRIC mapping remain outstanding.

## Validation

Unit tests verify module arithmetic, required gene coverage, finite/unique inputs, and the non-executable status of the proposed parameter links. The code is a feature-engineering interface and hypothesis registry only; it is not evidence of biological calibration or clinical utility.
