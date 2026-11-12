#!/usr/bin/env python3
"""Link I-SPY 2 MRI features to clinical outcomes and emit aggregate-only results.

The linked participant-level table is kept in memory and is never written.
Only aggregate JSON is written to the requested output path.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.impute import SimpleImputer


def exact_one_to_one_link(features: pd.DataFrame, clinical: pd.DataFrame) -> pd.DataFrame:
    """Perform a strict exact-ID one-to-one link; fail closed on ambiguous IDs."""
    f = features.copy()
    c = clinical.copy()
    fkey = "CLINICAL-TRIAL-SUBJECT-ID"
    ckey = "Patient_ID"
    required_f = {fkey, "VOLUME_TUM_BLU_V10", "VOLUME_TUM_BLU_V20",
                  "VOLUME_TUM_BLU_V30", "VOLUME_TUM_BLU_V40"}
    required_c = {ckey, "Arm", "HR", "HER2", "MP", "pCR", "Age_at_Screening"}
    if not required_f.issubset(f.columns):
        raise ValueError(f"Feature workbook missing columns: {sorted(required_f-set(f.columns))}")
    if not required_c.issubset(c.columns):
        raise ValueError(f"Clinical workbook missing columns: {sorted(required_c-set(c.columns))}")
    f[fkey] = pd.to_numeric(f[fkey], errors="raise").astype("Int64")
    c[ckey] = pd.to_numeric(c[ckey], errors="raise").astype("Int64")
    if f[fkey].isna().any() or c[ckey].isna().any():
        raise ValueError("Participant IDs cannot be missing")
    if f[fkey].duplicated().any() or c[ckey].duplicated().any():
        raise ValueError("Duplicate participant IDs: refusing ambiguous linkage")
    merged = f.merge(c, left_on=fkey, right_on=ckey, how="left",
                     validate="one_to_one", indicator=True)
    if not (merged["_merge"] == "both").all():
        missing = int((merged["_merge"] != "both").sum())
        raise ValueError(f"{missing} MRI participants lack clinical matches")
    merged = merged.drop(columns=["_merge"])
    return merged


def _or_row(y: pd.Series, x: pd.Series) -> dict:
    """Unadjusted OR per SD increase in x, via stable sklearn-free 2x1 GLM."""
    import statsmodels.api as sm
    x = pd.to_numeric(x, errors="raise").astype(float)
    z = (x - x.mean()) / x.std(ddof=0)
    fit = sm.Logit(y.astype(float), sm.add_constant(z)).fit(disp=False)
    b, se = float(fit.params.iloc[1]), float(fit.bse.iloc[1])
    return {"odds_ratio_per_1SD_increase": float(np.exp(b)),
            "ci95": [float(np.exp(b-1.96*se)), float(np.exp(b+1.96*se))],
            "p_value": float(fit.pvalues.iloc[1])}


def analyze(features: pd.DataFrame, clinical: pd.DataFrame) -> dict:
    d = exact_one_to_one_link(features, clinical)
    vcols = ["VOLUME_TUM_BLU_V10", "VOLUME_TUM_BLU_V20",
             "VOLUME_TUM_BLU_V30", "VOLUME_TUM_BLU_V40"]
    for col in vcols:
        d[col] = pd.to_numeric(d[col], errors="raise")
        if (d[col] <= 0).any():
            raise ValueError(f"Non-positive volume in {col}; log-ratio analysis invalid")
    for visit, col in zip(["T0", "T1", "T2", "T3"], vcols):
        d[f"logvol_{visit}"] = np.log(d[col])
    for visit in ["T1", "T2", "T3"]:
        d[f"logratio_{visit}_T0"] = d[f"logvol_{visit}"] - d["logvol_T0"]
        d[f"pctchange_{visit}"] = 100 * np.expm1(d[f"logratio_{visit}_T0"])

    y = pd.to_numeric(d["pCR"], errors="raise").astype(int)
    if set(y.unique()) != {0, 1}:
        raise ValueError("pCR must contain both binary codes 0 and 1")
    out = {
        "linkage": {
            "feature_n": int(len(features)),
            "clinical_n": int(len(clinical)),
            "matched_n": int(len(d)),
            "feature_unique_ids": int(features["CLINICAL-TRIAL-SUBJECT-ID"].nunique()),
            "clinical_unique_ids": int(clinical["Patient_ID"].nunique()),
            "feature_unmatched_n": 0,
            "clinical_only_n": int(clinical["Patient_ID"].nunique() -
                                   d["Patient_ID"].nunique()),
            "join": "exact one-to-one CLINICAL-TRIAL-SUBJECT-ID = Patient_ID"
        },
        "outcome": {
            "pCR_code_1_n": int((y == 1).sum()),
            "pCR_code_0_n": int((y == 0).sum()),
            "pCR_code_1_rate": float(y.mean()),
            "interpretation": "Code 1 is provisionally treated as pCR-positive; confirm against official data dictionary."
        },
        "pCR_by_HR_HER2": [],
        "pCR_by_arm": [],
        "ftv_candidate_change_by_pCR": {},
        "unadjusted_logistic": {},
        "internal_repeated_cv": {},
        "limitations": [
            "pCR code direction must be confirmed against the official data dictionary.",
            "No MRI acquisition dates; visits are not elapsed-day intervals.",
            "Treatment-arm comparisons are descriptive, not causal efficacy estimates.",
            "The workbook volume field's exact physical units and equivalence to trial FTV remain to be confirmed.",
            "Internal cross-validation is not external validation or clinical utility evidence.",
            "T3 is post-neoadjuvant/pre-surgery and is not an early prediction feature."
        ]
    }
    for keys, dest in [(["HR", "HER2"], "pCR_by_HR_HER2"), (["Arm"], "pCR_by_arm")]:
        grouped = d.groupby(keys, dropna=False)["pCR"].agg(["size", "sum", "mean"]).reset_index()
        out[dest] = [
            {**{k: (None if pd.isna(row[k]) else row[k].item() if hasattr(row[k], "item") else row[k])
                for k in keys},
             "n": int(row["size"]), "pCR_code_1_n": int(row["sum"]),
             "pCR_code_1_rate": float(row["mean"])}
            for _, row in grouped.iterrows()
        ]
    for code in [0, 1]:
        out["ftv_candidate_change_by_pCR"][str(code)] = {
            visit: {"median_pct_change": float(d.loc[y == code, f"pctchange_{visit}"].median()),
                    "mean_pct_change": float(d.loc[y == code, f"pctchange_{visit}"].mean())}
            for visit in ["T1", "T2", "T3"]
        }
    for visit in ["T1", "T2", "T3"]:
        a = d.loc[y == 1, f"logratio_{visit}_T0"]
        b = d.loc[y == 0, f"logratio_{visit}_T0"]
        out["unadjusted_logistic"][visit] = _or_row(y, d[f"logratio_{visit}_T0"])
        out["unadjusted_logistic"][visit]["mann_whitney_p"] = float(
            mannwhitneyu(a, b, alternative="two-sided").pvalue
        )

    # Repeated stratified CV is internal only. T1/T2 are separated from post-treatment T3.
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=20261112)
    for visit in ["T1", "T2"]:
        for with_arm in [False, True]:
            numeric = [f"logratio_{visit}_T0", "Age_at_Screening", "HR", "HER2"]
            transformers = [("num", Pipeline([
                ("impute", SimpleImputer(strategy="median")),
                ("scale", StandardScaler())
            ]), numeric)]
            cols = numeric + (["Arm"] if with_arm else [])
            if with_arm:
                transformers.append(("arm", OneHotEncoder(handle_unknown="ignore"), ["Arm"]))
            model = Pipeline([
                ("prep", ColumnTransformer(transformers)),
                ("clf", LogisticRegression(C=0.5, max_iter=3000))
            ])
            scores = cross_val_score(model, d[cols], y, cv=cv, scoring="roc_auc")
            key = f"{visit}_age_HR_HER2" + ("_plus_arm" if with_arm else "")
            out["internal_repeated_cv"][key] = {
                "mean_auc": float(scores.mean()),
                "fold_sd": float(scores.std()),
                "folds": int(len(scores))
            }
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", required=True)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--output", required=True,
                        help="Aggregate-only JSON output; no participant rows are written")
    args = parser.parse_args()
    features = pd.read_excel(args.features)
    clinical = pd.read_excel(args.clinical)
    result = analyze(features, clinical)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(result, indent=2, default=str) + "\n")
    print(json.dumps(result["linkage"], indent=2))


if __name__ == "__main__":
    main()
