# Data dictionary and integration contract

## Prototype model fields
| Field | Symbol | Meaning | Source |
|---|---|---|---|
| Cancer | C | Tumour-cell density proxy | Synthetic |
| Immune | I | Aggregate immune density proxy | Synthetic |
| Fibroblast | F | Stromal-cell proxy | Synthetic |
| Oxygen | O | Oxygen availability proxy | Synthetic |
| Nutrient | N | Nutrient availability proxy | Synthetic |
| Drug | D | Simulated compound field | Synthetic |
| Evasion | E | Immune-evasion activity proxy | Synthetic |
| Vasculature | V | Vascular density proxy | Synthetic |

## Future observation metadata
Each observation should carry cohort_id, pseudonymous patient_id, specimen_id, sample_id, assay_type, modality, platform, collection time, treatment context, endpoint definition where relevant, source accession, checksum, data-use permission and preprocessing version. Spatial data additionally need image geometry, pixel size, channel/marker map, segmentation method and mask provenance.

## Linkage rules
- Do not link samples to a patient without explicit source evidence.
- Do not treat unpaired cohort modalities as matched.
- Record missing and ambiguous values; never silently infer clinical outcomes.
- Keep endpoint and follow-up-window definitions explicit.
- Separate prognostic tasks from treatment-predictive tasks.
- Spatial fields require spatial data or must be labelled simulated.

## Protected data
Keep raw patient-level or access-controlled files in approved storage. Do not upload them to GitHub. Commit schemas, permitted metadata manifests and aggregate audits only.
