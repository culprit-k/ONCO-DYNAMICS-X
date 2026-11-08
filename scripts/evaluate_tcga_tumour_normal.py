"""Evaluate a leakage-safe TCGA-BRCA paired tumour/normal expression baseline."""
from __future__ import annotations

import argparse

from onco_dynamics_x.data_engine.tcga_baseline import evaluate_tumour_normal_baseline


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--top-k", type=int, default=100)
    args = parser.parse_args()
    report = evaluate_tumour_normal_baseline(
        args.matrix,
        args.manifest,
        args.output,
        random_seed=args.seed,
        n_splits=args.folds,
        top_k_genes=args.top_k,
    )
    print(f"Patients: {report['n_patients']}; samples: {report['n_samples']}")
    print(
        "Grouped holdout: balanced accuracy "
        f"{report['patient_grouped_holdout']['balanced_accuracy']:.3f}; "
        f"ROC-AUC {report['patient_grouped_holdout']['roc_auc']:.3f}"
    )
    print(
        "Grouped CV: balanced accuracy "
        f"{report['stratified_group_cv']['balanced_accuracy_mean']:.3f} ± "
        f"{report['stratified_group_cv']['balanced_accuracy_std']:.3f}; "
        f"ROC-AUC {report['stratified_group_cv']['roc_auc_mean']:.3f} ± "
        f"{report['stratified_group_cv']['roc_auc_std']:.3f}"
    )


if __name__ == "__main__":
    main()
