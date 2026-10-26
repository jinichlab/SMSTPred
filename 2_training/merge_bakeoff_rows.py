#!/usr/bin/env python3
"""Combine a bake-off re-run with the committed bake-off results.

The re-run (sam_classifier_comparison.py with BAKEOFF_ONLY=mlp) retrains only the MLP
variants. Because it replays the random stream, its train/test splits are identical to the committed run, so
the committed results of every other classifier (XGBoost, RF, LR, SVMs, MMseqs2, Foldseek, Profile-HMM) still
apply. This script writes a COPY of the committed detailed / per-label / summary CSVs in which the MLP rows of
the chosen dataset(s) are replaced by the re-run's; everything else is untouched. --keep-datasets drops the other
datasets' rows from the output (used to keep only the full-set results).

  python merge_bakeoff_rows.py --base <original bake-off folder> --new <MLP re-run folder> \
                               --out <merged folder> --datasets full --keep-datasets full

--classifiers limits the replacement to the named MLP rows, and --base may be a previous merge (an absolute path
works too), so single variants can be re-run and swapped in afterwards. See
sam_classifier_comparison_output/PROVENANCE.md for how the committed results were assembled.
"""
import argparse
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
METRICS = ["auprc", "auroc", "mcc", "f1", "accuracy"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="sam_classifier_comparison_output", help="committed results folder")
    ap.add_argument("--new", required=True, help="re-run folder (MLP rows only)")
    ap.add_argument("--out", required=True, help="folder to write the merged copy to")
    ap.add_argument("--datasets", nargs="+", default=["full"], help="datasets whose MLP rows are replaced")
    ap.add_argument("--classifiers", nargs="+", default=None,
                    help="replace only these classifiers' rows (default: every MLP_ variant in the re-run)")
    ap.add_argument("--keep-datasets", nargs="+", default=None, help="write only these datasets' rows (default: all)")
    a = ap.parse_args()
    base, new, out = (HERE / a.base, HERE / a.new, HERE / a.out)
    out.mkdir(parents=True, exist_ok=True)

    for name in ("detailed", "per_label"):
        f = f"sam_classifier_comparison_{name}.csv"
        B, N = pd.read_csv(base / f), pd.read_csv(new / f)
        chosen = (lambda df: df.classifier.isin(a.classifiers)) if a.classifiers else (lambda df: df.classifier.str.startswith("MLP_"))
        drop = chosen(B) & B.dataset.isin(a.datasets)
        take = chosen(N) & N.dataset.isin(a.datasets)
        merged = pd.concat([B[~drop], N[take]], ignore_index=True)
        merged = merged.sort_values(["split_mode", "dataset", "iteration", "classifier"]).reset_index(drop=True)
        assert len(merged) == len(B), (f, len(merged), len(B))          # same rows, only values replaced
        if a.keep_datasets:
            merged = merged[merged.dataset.isin(a.keep_datasets)].reset_index(drop=True)
        merged.to_csv(out / f, index=False)
        print(f"{f}: {int(drop.sum())} committed MLP rows replaced by {int(take.sum())} re-run rows -> {len(merged)} rows")

    det = pd.read_csv(out / "sam_classifier_comparison_detailed.csv")
    summ = det.groupby(["split_mode", "dataset", "classifier"])[METRICS].agg(["mean", "std"]).reset_index()
    summ.columns = ["_".join(c).strip("_") for c in summ.columns]
    summ.to_csv(out / "sam_classifier_comparison_summary.csv", index=False)
    print(f"summary recomputed from the merged detailed table -> {out}")


if __name__ == "__main__":
    main()
