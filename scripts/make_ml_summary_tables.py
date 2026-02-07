#!/usr/bin/env python3
"""make_ml_summary_tables.py

Generate two tables from machine learning benchmark results:

1) Main Table X: concise, reviewer-friendly summary.
2) Supplementary Table Sx: full ranked benchmark (e.g., 113 combinations).

Inputs are CSVs. Column names are normalized/aliased to be robust to naming
inconsistencies (e.g., train_cv_auc/cv_auc/auc_train_cv, etc.).

Outputs are written into a new folder (default: ./ml_summary_tables):
- full_benchmark_ranked.csv
- full_benchmark_ranked.xlsx
- main_table_x.csv
- main_table_x.xlsx
- main_table_x.md

Example:
  python scripts/make_ml_summary_tables.py \
    --benchmark_csv results/benchmark_113.csv \
    --lmbd_csv results/lmbd_final.csv \
    --single_pathway_csv results/single_pathway.csv \
    --out_dir ml_summary_tables
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


def _require_pandas() -> "tuple[object, object]":
    try:
        import numpy as np  # type: ignore
        import pandas as pd  # type: ignore
    except Exception as e:
        raise SystemExit(
            "This script requires pandas and numpy. "
            "Install them (e.g., pip install pandas numpy openpyxl).\n"
            f"Import error: {e}"
        )
    return pd, np


def _norm_col(s: str) -> str:
    s = str(s).strip().lower()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def _first_present_col(df, aliases: Sequence[str]) -> Optional[str]:
    """Return the first matching column in df, given alias candidates.

    Matching is done on normalized column names.
    """
    norm_to_real: Dict[str, str] = {_norm_col(c): c for c in df.columns}
    for a in aliases:
        c = norm_to_real.get(_norm_col(a))
        if c is not None:
            return c
    return None


def _to_num(series):
    pd, _np = _require_pandas()
    return pd.to_numeric(series, errors="coerce")


def _safe_div(num, den):
    pd, np = _require_pandas()
    num = _to_num(num)
    den = _to_num(den)
    out = num / den
    out = out.where(den != 0)
    return out


def _canonical_algorithm(val: object) -> str:
    """Map raw algorithm text to a canonical label for grouping."""
    s = "" if val is None else str(val)
    s0 = s
    s = s.strip().lower()
    s = re.sub(r"\s+", " ", s)

    # Common patterns (keep conservative; fall back to original string).
    patterns: List[Tuple[str, str]] = [
        (r"\b(random\s*forest|rand\s*forest|rf)\b", "RF"),
        (r"\b(xgboost|xgb)\b", "XGBoost"),
        (r"\b(lightgbm|lgbm)\b", "GBM"),
        (r"\b(gradient\s*boost|gbm|gboost)\b", "GBM"),
        (r"\b(svm|support\s*vector)\b", "SVM"),
        (r"\b(naive\s*bayes|naivebayes|\bnb\b)\b", "NB"),
        (r"\b(pls|partial\s*least\s*squares)\b", "PLS"),
        (r"\b(dnn|mlp|neural\s*net|deep\s*learning)\b", "DNN"),
        (r"\b(elastic\s*net|elasticnet|enet)\b", "ElasticNet"),
        (r"\b(lasso)\b", "Lasso"),
        (r"\b(ridge)\b", "Ridge"),
        (r"\b(logistic\s*regression|logit|\blr\b)\b", "LR"),
        (r"\b(gmm\s*\-?\s*lr|gmm\s*lr)\b", "GMM-LR"),
    ]
    for pat, canon in patterns:
        if re.search(pat, s, flags=re.IGNORECASE):
            return canon

    # Heuristic: if model string contains algo keywords.
    if re.search(r"\b(random\s*forest|\brf\b)", s, flags=re.IGNORECASE):
        return "RF"
    if re.search(r"\b(xgboost|\bxgb\b)", s, flags=re.IGNORECASE):
        return "XGBoost"
    if re.search(r"\bsvm\b", s, flags=re.IGNORECASE):
        return "SVM"

    # Fallback.
    s0 = str(s0).strip()
    return s0 if s0 else "Unknown"


def _read_csv_flexible(path: Path):
    pd, _np = _require_pandas()
    # Try utf-8-sig first (handles BOM), then fall back to default.
    try:
        return pd.read_csv(path, encoding="utf-8-sig")
    except UnicodeDecodeError:
        return pd.read_csv(path)


def _standardize_ml_table(df, *, source_label: str) -> "object":
    """Standardize columns for benchmark/lmbd/single-pathway tables.

    Returns a new DataFrame with standardized columns:
      Model, Algorithm, FeatureSet,
      Train_CV_AUC, External_AUC, Mean_AUC,
      Accuracy, Sensitivity, Specificity,
      Notes

    Also includes helper columns:
      _src_index, _algo_canon
    """
    pd, np = _require_pandas()
    df0 = df.copy()
    df0["_src_index"] = range(len(df0))

    # Core identity columns.
    model_col = _first_present_col(
        df0,
        [
            "Model",
            "ModelName",
            "model_name",
            "name",
            "pipeline",
            "estimator",
            "classifier",
        ],
    )
    algo_col = _first_present_col(
        df0,
        [
            "Algorithm",
            "algo",
            "method",
            "model_type",
            "classifier",
        ],
    )
    feat_col = _first_present_col(
        df0,
        [
            "FeatureSet",
            "feature_set",
            "featureset",
            "features",
            "feature",
            "signature",
            "pathway",
            "pathway_name",
            "pathway_model",
        ],
    )

    # AUC columns (aliases).
    train_auc_col = _first_present_col(
        df0,
        [
            "Train_CV_AUC",
            "train_cv_auc",
            "cv_auc",
            "auc_train_cv",
            "auc_cv",
            "train_auc",
            "auc_train",
        ],
    )
    ext_auc_col = _first_present_col(
        df0,
        [
            "External_AUC",
            "external_auc",
            "test_auc",
            "val_auc",
            "auc_external",
            "auc_test",
            "auc_val",
        ],
    )
    mean_auc_col = _first_present_col(
        df0,
        [
            "Mean_AUC",
            "mean_auc",
            "avg_auc",
            "auc_mean",
            "meanAUC",
        ],
    )

    # Other metrics.
    acc_col = _first_present_col(df0, ["Accuracy", "accuracy", "acc"])
    sens_col = _first_present_col(
        df0, ["Sensitivity", "sensitivity", "recall", "tpr", "true_positive_rate"]
    )
    spec_col = _first_present_col(
        df0, ["Specificity", "specificity", "tnr", "true_negative_rate"]
    )
    tp_col = _first_present_col(df0, ["tp", "TP", "true_positive"])
    tn_col = _first_present_col(df0, ["tn", "TN", "true_negative"])
    fp_col = _first_present_col(df0, ["fp", "FP", "false_positive"])
    fn_col = _first_present_col(df0, ["fn", "FN", "false_negative"])

    out = pd.DataFrame()
    out["Model"] = df0[model_col] if model_col else ""

    if algo_col:
        out["Algorithm"] = df0[algo_col]
    else:
        # Try to infer algorithm from model column.
        out["Algorithm"] = out["Model"].map(_canonical_algorithm)

    out["FeatureSet"] = df0[feat_col] if feat_col else ""

    out["Train_CV_AUC"] = _to_num(df0[train_auc_col]) if train_auc_col else np.nan
    out["External_AUC"] = _to_num(df0[ext_auc_col]) if ext_auc_col else np.nan
    out["Mean_AUC"] = _to_num(df0[mean_auc_col]) if mean_auc_col else np.nan

    # Compute Mean_AUC if missing.
    mean_missing = out["Mean_AUC"].isna()
    can_compute = (~out["Train_CV_AUC"].isna()) & (~out["External_AUC"].isna())
    out.loc[mean_missing & can_compute, "Mean_AUC"] = (
        out.loc[mean_missing & can_compute, "Train_CV_AUC"]
        + out.loc[mean_missing & can_compute, "External_AUC"]
    ) / 2.0

    out["Accuracy"] = _to_num(df0[acc_col]) if acc_col else np.nan
    out["Sensitivity"] = _to_num(df0[sens_col]) if sens_col else np.nan
    out["Specificity"] = _to_num(df0[spec_col]) if spec_col else np.nan

    # Compute metrics from confusion matrix if present and missing.
    if tp_col and tn_col and fp_col and fn_col:
        tp = _to_num(df0[tp_col])
        tn = _to_num(df0[tn_col])
        fp = _to_num(df0[fp_col])
        fn = _to_num(df0[fn_col])
        total = tp + tn + fp + fn
        if out["Accuracy"].isna().any():
            out.loc[out["Accuracy"].isna(), "Accuracy"] = _safe_div(tp + tn, total)
        if out["Sensitivity"].isna().any():
            out.loc[out["Sensitivity"].isna(), "Sensitivity"] = _safe_div(tp, tp + fn)
        if out["Specificity"].isna().any():
            out.loc[out["Specificity"].isna(), "Specificity"] = _safe_div(tn, tn + fp)

    out["Notes"] = ""
    out["_src_index"] = df0["_src_index"]
    out["_algo_canon"] = out["Algorithm"].map(_canonical_algorithm)

    # Attach source label for easier debugging.
    out["_source"] = source_label

    return out


def _rank_benchmark(df_std):
    pd, np = _require_pandas()
    df = df_std.copy()

    def _sort_key(s):
        s2 = _to_num(s)
        return s2.fillna(-np.inf)

    df = df.sort_values(
        by=["Mean_AUC", "External_AUC", "Train_CV_AUC"],
        ascending=[False, False, False],
        key=_sort_key,
        kind="mergesort",  # stable
    ).reset_index(drop=True)
    df.insert(0, "Rank", range(1, len(df) + 1))
    return df


def _select_winner_row(df_ranked):
    """Select benchmark winner; if tied, prefer RF among ties."""
    pd, np = _require_pandas()
    if len(df_ranked) == 0:
        return None
    best = _to_num(df_ranked["Mean_AUC"]).max(skipna=True)
    if pd.isna(best):
        return df_ranked.iloc[0]
    ties = df_ranked[np.isclose(_to_num(df_ranked["Mean_AUC"]), best, equal_nan=False)]
    if len(ties) == 0:
        return df_ranked.iloc[0]
    rf_ties = ties[ties["_algo_canon"].astype(str).str.upper().eq("RF")]
    if len(rf_ties) > 0:
        return rf_ties.iloc[0]
    return ties.iloc[0]


def _pick_best_row(df_std, *, prefer: Sequence[str] = ("External_AUC", "Mean_AUC", "Train_CV_AUC")):
    pd, np = _require_pandas()
    if len(df_std) == 0:
        return None

    df = df_std.copy()
    keys = list(prefer)
    # Sort with NaNs at bottom.
    for k in keys:
        if k not in df.columns:
            df[k] = np.nan
    df = df.sort_values(
        by=keys,
        ascending=[False] * len(keys),
        key=lambda s: _to_num(s).fillna(-np.inf),
        kind="mergesort",
    )
    return df.iloc[0]


def _build_main_table(
    benchmark_ranked,
    *,
    lmbd_std=None,
    single_pathway_std=None,
    strong_algo_set: Sequence[str],
    strong_min: int,
    strong_max: int,
    single_pathway_max: int,
):
    pd, np = _require_pandas()

    cols = [
        "Group",
        "Model",
        "Algorithm",
        "FeatureSet",
        "Train_CV_AUC",
        "External_AUC",
        "Mean_AUC",
        "Accuracy",
        "Sensitivity",
        "Specificity",
        "Notes",
    ]

    picked_src: set[Tuple[str, int]] = set()
    rows: List[Dict[str, object]] = []

    def _add_from_row(group: str, row, notes: str):
        if row is None:
            return
        src = str(row.get("_source", ""))
        idx = int(row.get("_src_index", -1))
        key = (src, idx)
        if key in picked_src:
            return
        picked_src.add(key)
        rows.append(
            {
                "Group": group,
                "Model": row.get("Model", ""),
                "Algorithm": row.get("Algorithm", ""),
                "FeatureSet": row.get("FeatureSet", ""),
                "Train_CV_AUC": row.get("Train_CV_AUC", np.nan),
                "External_AUC": row.get("External_AUC", np.nan),
                "Mean_AUC": row.get("Mean_AUC", np.nan),
                "Accuracy": row.get("Accuracy", np.nan),
                "Sensitivity": row.get("Sensitivity", np.nan),
                "Specificity": row.get("Specificity", np.nan),
                "Notes": notes,
            }
        )

    # (1) Benchmark winner.
    winner = _select_winner_row(benchmark_ranked)
    _add_from_row(
        "Benchmark winner",
        winner,
        "best in 113 benchmark by mean AUC",
    )

    # (2) Final model (LMBD, GMM-LR).
    if lmbd_std is not None and len(lmbd_std) > 0:
        lmbd_best = _pick_best_row(lmbd_std, prefer=("External_AUC", "Mean_AUC", "Train_CV_AUC"))
        if lmbd_best is not None:
            # Ensure Algorithm is present and recognizable.
            lmbd_best = lmbd_best.copy()
            algo_c = _canonical_algorithm(lmbd_best.get("Algorithm", ""))
            if algo_c == "Unknown":
                lmbd_best["Algorithm"] = "GMM-LR"
            _add_from_row(
                "Final model (LMBD, GMM-LR)",
                lmbd_best,
                "final interpretable model",
            )

    # (3) Strong baselines (3-6): best per algorithm among benchmark.
    if len(benchmark_ranked) > 0:
        bench = benchmark_ranked.copy()
        bench = bench[~bench[["_source", "_src_index"]].apply(tuple, axis=1).isin(picked_src)]
        bench["_algo_canon"] = bench["Algorithm"].map(_canonical_algorithm)
        bench = bench[bench["_algo_canon"].isin(list(strong_algo_set))]

        best_per_algo = []
        for algo in strong_algo_set:
            sub = bench[bench["_algo_canon"].eq(algo)]
            if len(sub) == 0:
                continue
            best_row = _pick_best_row(sub, prefer=("Mean_AUC", "External_AUC", "Train_CV_AUC"))
            if best_row is not None:
                best_per_algo.append(best_row)
        if best_per_algo:
            best_df = pd.DataFrame(best_per_algo)
            best_df = best_df.sort_values(
                by=["Mean_AUC", "External_AUC", "Train_CV_AUC"],
                ascending=[False, False, False],
                key=lambda s: _to_num(s).fillna(-np.inf),
                kind="mergesort",
            )
            best_df = best_df.head(strong_max)
            for _, r in best_df.iterrows():
                _add_from_row(
                    "Strong baselines",
                    r,
                    f"best {r.get('_algo_canon', r.get('Algorithm', ''))} baseline (by mean AUC)",
                )

    # (4) Single-pathway baselines (top by external_auc).
    if single_pathway_std is not None and len(single_pathway_std) > 0:
        sp = single_pathway_std.copy()
        sp_best = sp.sort_values(
            by=["External_AUC", "Mean_AUC", "Train_CV_AUC"],
            ascending=[False, False, False],
            key=lambda s: _to_num(s).fillna(-np.inf),
            kind="mergesort",
        ).head(single_pathway_max)
        for _, r in sp_best.iterrows():
            _add_from_row(
                "Single-pathway baselines",
                r,
                "single-pathway baseline (ranked by external AUC)",
            )

    main = pd.DataFrame(rows, columns=cols)

    # Enforce 3-6 strong baseline rows if possible; do not drop other groups.
    if "Strong baselines" in set(main["Group"].astype(str)):
        strong_rows = main[main["Group"].eq("Strong baselines")]
        if len(strong_rows) < strong_min:
            # Keep as-is; caller may have limited algorithms present.
            pass

    return main


def _df_to_markdown(df, *, float_digits: int = 3) -> str:
    # Lightweight markdown generation without extra dependencies.
    def fmt(v):
        if v is None:
            return ""
        try:
            import math

            if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
                return ""
        except Exception:
            pass
        # Avoid printing pandas NaN.
        try:
            pd, _np = _require_pandas()
            if pd.isna(v):
                return ""
        except Exception:
            pass
        if isinstance(v, (int,)):
            return str(v)
        if isinstance(v, (float,)):
            return f"{v:.{float_digits}f}"
        return str(v)

    cols = list(df.columns)
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join(["---"] * len(cols)) + " |"
    lines = [header, sep]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(fmt(r[c]) for c in cols) + " |")
    return "\n".join(lines) + "\n"


def _write_outputs(df, *, out_dir: Path, csv_name: str, xlsx_name: str, sheet_name: str):
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / csv_name
    xlsx_path = out_dir / xlsx_name
    df.to_csv(csv_path, index=False)
    try:
        df.to_excel(xlsx_path, index=False, sheet_name=sheet_name)
    except Exception as e:
        print(
            f"WARNING: Failed to write Excel file {xlsx_path}: {e}\n"
            "Install openpyxl to enable .xlsx output.",
            file=sys.stderr,
        )


def main(argv: Optional[Sequence[str]] = None) -> int:
    pd, np = _require_pandas()

    p = argparse.ArgumentParser(
        description="Generate ML summary tables (Table X + Supplementary Table Sx) from benchmark CSVs.",
    )
    p.add_argument("--benchmark_csv", required=True, type=Path, help="Benchmark CSV (e.g., 113 combinations).")
    p.add_argument("--lmbd_csv", type=Path, default=None, help="Optional LMBD (GMM-LR) CSV.")
    p.add_argument("--single_pathway_csv", type=Path, default=None, help="Optional single-pathway baseline CSV.")
    p.add_argument(
        "--out_dir",
        type=Path,
        default=Path("ml_summary_tables"),
        help="Output directory (created if missing).",
    )
    p.add_argument("--sheet_name", default="Table", help="Excel sheet name.")
    p.add_argument("--float_digits", type=int, default=3, help="Digits for markdown float formatting.")

    p.add_argument(
        "--lmbd_id_col",
        default=None,
        help="Optional ID column in LMBD CSV to select a specific row.",
    )
    p.add_argument(
        "--lmbd_id_value",
        default=None,
        help="Optional ID value in LMBD CSV to select a specific row.",
    )

    p.add_argument(
        "--strong_algorithms",
        default="LR,Ridge,Lasso,ElasticNet,SVM,GBM,XGBoost,NB,PLS,DNN",
        help="Comma-separated canonical algorithm set for strong baselines.",
    )
    p.add_argument("--strong_min", type=int, default=3, help="Minimum strong baseline rows (best-effort).")
    p.add_argument("--strong_max", type=int, default=6, help="Maximum strong baseline rows.")
    p.add_argument(
        "--single_pathway_topk",
        type=int,
        default=3,
        help="Top-K single-pathway baselines by external AUC.",
    )

    args = p.parse_args(list(argv) if argv is not None else None)

    if not args.benchmark_csv.exists():
        raise SystemExit(f"benchmark_csv not found: {args.benchmark_csv}")

    bench_raw = _read_csv_flexible(args.benchmark_csv)
    bench_std = _standardize_ml_table(bench_raw, source_label="benchmark")

    # Basic validation: require either mean AUC or both train/external AUC.
    if bench_std["Mean_AUC"].isna().all() and (
        bench_std["Train_CV_AUC"].isna().all() or bench_std["External_AUC"].isna().all()
    ):
        raise SystemExit(
            "Could not find usable AUC columns in benchmark_csv. "
            "Need mean_auc OR both (train_cv_auc and external_auc) via supported aliases."
        )

    bench_ranked = _rank_benchmark(bench_std)

    # Supplementary Table Sx.
    full_cols = [
        "Rank",
        "Model",
        "Algorithm",
        "FeatureSet",
        "Train_CV_AUC",
        "External_AUC",
        "Mean_AUC",
        "Accuracy",
        "Sensitivity",
        "Specificity",
        "Notes",
    ]
    full_out = bench_ranked[full_cols].copy()
    _write_outputs(
        full_out,
        out_dir=args.out_dir,
        csv_name="full_benchmark_ranked.csv",
        xlsx_name="full_benchmark_ranked.xlsx",
        sheet_name=args.sheet_name,
    )

    # Optional LMBD.
    lmbd_std = None
    if args.lmbd_csv is not None:
        if not args.lmbd_csv.exists():
            raise SystemExit(f"lmbd_csv not found: {args.lmbd_csv}")
        lmbd_raw = _read_csv_flexible(args.lmbd_csv)
        lmbd_std = _standardize_ml_table(lmbd_raw, source_label="lmbd")
        # Optional selection by ID.
        if args.lmbd_id_col and args.lmbd_id_value:
            id_col = _first_present_col(lmbd_raw, [args.lmbd_id_col])
            if not id_col:
                raise SystemExit(f"lmbd_id_col not found in lmbd_csv: {args.lmbd_id_col}")
            sel = lmbd_raw[id_col].astype(str) == str(args.lmbd_id_value)
            lmbd_std = lmbd_std.loc[sel.values].copy()
            if len(lmbd_std) == 0:
                raise SystemExit(
                    f"No LMBD rows matched {args.lmbd_id_col}={args.lmbd_id_value}"
                )

        # Force recognizable algorithm label for display.
        lmbd_std.loc[:, "Algorithm"] = lmbd_std["Algorithm"].apply(
            lambda x: "GMM-LR" if _canonical_algorithm(x) in ("Unknown", "") else x
        )
        lmbd_std.loc[:, "_algo_canon"] = lmbd_std["Algorithm"].map(_canonical_algorithm)

    # Optional single-pathway.
    sp_std = None
    if args.single_pathway_csv is not None:
        if not args.single_pathway_csv.exists():
            raise SystemExit(f"single_pathway_csv not found: {args.single_pathway_csv}")
        sp_raw = _read_csv_flexible(args.single_pathway_csv)
        sp_std = _standardize_ml_table(sp_raw, source_label="single_pathway")
        # If algorithm missing, label it.
        sp_std.loc[:, "Algorithm"] = sp_std["Algorithm"].replace({"": "Single-pathway"})

    strong_algo_set = [s.strip() for s in str(args.strong_algorithms).split(",") if s.strip()]

    main_df = _build_main_table(
        bench_ranked,
        lmbd_std=lmbd_std,
        single_pathway_std=sp_std,
        strong_algo_set=strong_algo_set,
        strong_min=args.strong_min,
        strong_max=args.strong_max,
        single_pathway_max=args.single_pathway_topk,
    )

    _write_outputs(
        main_df,
        out_dir=args.out_dir,
        csv_name="main_table_x.csv",
        xlsx_name="main_table_x.xlsx",
        sheet_name=args.sheet_name,
    )

    md_path = args.out_dir / "main_table_x.md"
    md_path.write_text(_df_to_markdown(main_df, float_digits=args.float_digits), encoding="utf-8")

    print(f"Wrote outputs to: {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
