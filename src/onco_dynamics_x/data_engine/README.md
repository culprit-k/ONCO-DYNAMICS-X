# Data engine

Implemented:
- `manifest.py`: required metadata columns, unique sample IDs, access-label checks and warnings for controlled/restricted data.
- `expression.py`: identifier-aware samples-by-genes validation, variance filtering and train-only scaling with explicit gene reordering.

Planned: cohort downloads, assay-specific QC, modality pairing, image/mask metadata and approved controlled-data workflows. Never commit restricted raw source data. Fit preprocessing on training data only and apply the frozen transform to validation/test cohorts.
