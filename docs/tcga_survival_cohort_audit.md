# TCGA-BRCA PanCancer Atlas survival-cohort audit

This workflow audits the cBioPortal study brca_tcga_pan_can_atlas_2018 as a candidate development cohort for an endpoint-specific survival model.

It collects aggregate study, patient/sample, molecular-profile and clinical-attribute metadata. The audit does not download expression matrices, export row-level clinical records, fit a model, or claim external validation.

Run locally:

    pip install -e ".[dev]"
    pytest tests/test_cbioportal_audit.py
    python scripts/audit_cbioportal_tcga.py

Or run the GitHub Actions workflow named TCGA-BRCA survival development cohort audit. Its artifact contains an aggregate JSON report and a clinical-attribute coverage CSV only.

The current paired TCGA tumour/normal benchmark is a separate 48-patient engineering experiment and is not the intended survival-model development cohort. The PanCancer Atlas audit must establish the available patient-level OS fields, event/time completeness, molecular profile assay and sample linkage before a training pipeline is designed.

The molecular assay scale must be read from the returned profile metadata. Do not assume compatibility with METABRIC Illumina HT-12 v3 microarray data. Any cross-cohort validation requires frozen preprocessing, common gene mapping and cohort/platform diagnostics.

No treatment-response inference is made from survival fields.
