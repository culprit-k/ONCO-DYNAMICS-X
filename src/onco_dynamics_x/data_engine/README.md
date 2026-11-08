# Data engine

Implemented:
- GDC public/open file metadata search, sample linkage export and selected open-file download with checksum verification.
- Strict expression table parsing and matrix assembly from downloaded files. Sample IDs remain attached to rows and gene columns are aligned by stable gene ID.
- Cohort manifest checks and train-only expression preprocessing.
- Patient-grouped train/test splitting.

The matrix assembler defaults to GDC STAR Counts column `tpm_unstranded`; verify the actual file headers and intended quantification scale for the selected workflow. It rejects mismatched gene sets and multi-sample linkage within a file rather than silently guessing.

Controlled, restricted or identifiable participant data must remain in approved storage and must not be committed to this public repository.
