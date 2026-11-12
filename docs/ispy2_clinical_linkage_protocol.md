# I-SPY 2 clinical metadata linkage protocol

The TCIA release contains at least two distinct resources relevant here: (1) the 384-participant multi-feature MRI workbook and (2) a separate clinical workbook/CSV for the 985-participant I-SPY 2 Imaging Cohort 1 (719 I-SPY 2 plus 266 ACRIN-6698 participants). The latter is not a verified row-level companion to the 384-participant subset. Do not infer shared row order, identical identifiers, or matching cohort membership.

## Required steps before any join

1. Download the clinical release from the current [TCIA I-SPY 2 collection page](https://www.cancerimagingarchive.net/collection/ispy2/) and record release/version, date, license and source citation.
2. Inspect the workbook/data dictionary to establish the documented linkage key. Do not join on names or row order.
3. Compare key type and formatting without exposing key values in reports. Check key uniqueness and missingness independently in both files.
4. Default to a one-to-one join. Select many-to-one only if the source documentation explicitly supports multiple feature rows per clinical subject and one clinical row per subject. Reject many-to-many joins.
5. Report feature rows, clinical rows, matched feature rows, unmatched feature rows, unmatched clinical rows and duplicate-key counts. Do not print or export identifiers in the audit report.
6. Inspect the clinical columns for treatment arm/regimen, regimen start/end dates, HR/ER/PR/HER2 definitions, pCR definition, MRI scan dates and the exact date/time origin. Verify units, specimen timing and field semantics.
7. Use actual acquisition dates to derive elapsed days relative to a documented reference. Nominal visit labels alone cannot support per-day kinetics.
8. Keep row-level linked data and all linkage tables in approved private storage. Commit only code, aggregate summaries and documentation to the public repository.

## Local linkage validator

After confirming the documented key names and join cardinality, run:

```bash
python scripts/validate_ispy2_linkage.py \
  --features /path/to/Multi-feature-MRI-NACT-Data.xlsx \
  --clinical /path/to/ISPY2-Clinical-Data.xlsx \
  --feature-key DOCUMENTED_FEATURE_KEY \
  --clinical-key DOCUMENTED_CLINICAL_KEY \
  --cardinality one_to_one \
  --report /path/to/private/ispy2_linkage_audit.json
```

The report contains counts, column names and overlap status but no key values or row-level records. The default mode does not write a merged table. If a local merged table is required, add `--linked-output /path/to/private/linked.csv`; keep that file local and do not commit it. Many-to-many joins are prohibited and duplicated clinical keys are rejected.

## Analysis gate

Do not fit treatment-adjusted models until the linkage is documented, duplicate and unmatched cases are explained, pCR/subtype/treatment fields are source-validated, and visit-specific acquisition times are available. A regimen field alone does not provide drug concentration or dosing history, and cannot identify a pharmacodynamic EC50.
