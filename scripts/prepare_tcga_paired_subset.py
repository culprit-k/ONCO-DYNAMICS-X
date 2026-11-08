"""Select paired primary tumour and normal TCGA samples from GDC metadata."""
from __future__ import annotations

import argparse

from onco_dynamics_x.data_engine.tcga_baseline import select_paired_samples


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tumour-metadata", required=True)
    parser.add_argument("--normal-metadata", required=True)
    parser.add_argument("--output-metadata", required=True)
    parser.add_argument("--output-manifest", required=True)
    parser.add_argument("--max-pairs", type=int, default=20)
    args = parser.parse_args()
    rows = select_paired_samples(
        args.tumour_metadata,
        args.normal_metadata,
        args.output_metadata,
        args.output_manifest,
        max_pairs=args.max_pairs,
    )
    print(
        f"Selected {len(rows) // 2} paired patients and {len(rows)} specimens "
        f"for patient-grouped evaluation."
    )


if __name__ == "__main__":
    main()
