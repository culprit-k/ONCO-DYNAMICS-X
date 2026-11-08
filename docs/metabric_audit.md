# METABRIC independent-cohort audit

## Objective

This audit checks whether the public cBioPortal METABRIC study has the cohort structure, molecular profiles and patient-level clinical fields needed to plan a defensible external prognostic validation.

The audit uses the public cBioPortal REST API v3 and exports aggregate counts plus a clinical-attribute coverage inventory. It does not export patient IDs, sample IDs, or individual clinical values.

## Reproduce

Install the development dependencies, run the focused tests, then execute:

    pip install -e ".[dev]"
    pytest tests/test_cbioportal_audit.py
    python scripts/audit_cbioportal_metabric.py

Or run the GitHub Actions workflow named METABRIC independent cohort audit. Its artifact contains:
- report.json: study metadata, aggregate cohort/linkage counts, molecular profile metadata, sample-list metadata, endpoint-field candidates and audit checks.
- clinical_attribute_coverage.csv: attribute identifiers/names, patient/sample level and aggregate non-missing coverage counts.

The API base URL defaults to https://www.cbioportal.org/api/v3. The client uses standard-library HTTP/JSON functionality and does not require credentials for public study data.

## What the audit checks

- Study metadata and cBioPortal study identifier.
- Unique patient and sample counts.
- Sample-to-patient linkage completeness.
- Molecular profile inventory, including profile names and data-type metadata.
- Sample-list inventory.
- Patient- and sample-level clinical attribute dictionary.
- Aggregate non-empty coverage for each clinical attribute.
- Name-based candidates for overall survival, disease-specific survival, relapse-free/disease-free survival, progression-free survival and treatment-related fields.

## Interpretation safeguards

A field whose name resembles OS_MONTHS or OS_STATUS is only a candidate. The script does not establish the endpoint's time origin, unit, event definition, censoring convention, follow-up completeness or clinical meaning. These must be checked against the study's data dictionary and source publication.

Likewise, a treatment-related attribute is not necessarily a treatment-response label. Treatment exposure, line/context, response criteria and assessment timing must be linked and verified before any response model is considered.

The molecular profile inventory is metadata only. A profile being called mRNA expression does not prove that it is compatible with TCGA STAR-count/TPM data. METABRIC assay platform, probe annotation, gene mapping and scale must be reviewed before model transfer.

The workflow deliberately does not download the full molecular matrix or fit a survival model. Those steps require a successful audit and a prespecified endpoint/harmonization plan.

## Data governance

Only aggregate audit outputs are uploaded as a short-lived workflow artifact. Raw clinical records and identifiers remain in the workflow process and are not committed. Do not add row-level clinical values or participant identifiers to this public repository.

## API references

- cBioPortal API and clients: https://docs.cbioportal.org/web-api-and-clients/
- cBioPortal API documentation: https://www.cbioportal.org/api
- Study: https://www.cbioportal.org/study/summary?id=brca_metabric
