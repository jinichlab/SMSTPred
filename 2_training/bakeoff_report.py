#!/usr/bin/env python3
"""
Stage 2 -- everything that is made from the 16 bake-off iterations (after run_bakeoff.sh):

  python bakeoff_report.py [--out <release>/benchmark]

  1. merges iterations/iterNN_*.csv into sam_classifier_comparison_{detailed,per_label,summary,per_label_summary}.csv (the files the figure / table / statistics scripts read)
  2. classifier_mcc_stats.py (Friedman, Nemenyi, Bonferroni-Dunn vs MLP_ProstT5, Wilcoxon-Holm), plot_cluster_mcc_zoom.py (Fig 2A), build_metrics_table.py (Table 1), run on this folder
  3. paired_tests_vs_baselines.csv: every MLP minus each search baseline, paired over the splits (mean difference, 95 % bootstrap interval, splits won, Wilcoxon p, Holm over that MLP's
     three comparisons). The 16 splits are random re-draws of the same proteins, so these p-values are optimistic (the splits overlap); read the differences and the split counts first.
  4. train_val_loss_by_embedding.png / train_val_loss_by_label.png and kept_epochs.csv: the train / validation loss over all 2,000 epochs and the epoch early stopping kept in every split
Prints the main tables. Needs numpy, pandas, scipy, matplotlib and scikit-posthocs (the conda environment).
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import wilcoxon

HERE = Path(__file__).resolve().parent
LABELS = ["Metabolite", "Protein", "RNA", "DNA"]
METRICS = ["auprc", "auroc", "mcc", "f1", "accuracy", "mcc_valtuned"]
EMBS = ["ProstT5", "SaProt", "ProtT5", "ESM-IF", "ProtT5+ESM-IF", "ProstT5+ESM-IF"]
SEARCH = ["MMseqs2_1NN", "Foldseek_1NN", "Profile_HMM"]


def holm(p):
    p = np.asarray(p, float)
    order, m, out, run = np.argsort(p), len(p), np.empty(len(p)), 0.0
    for r, i in enumerate(order):
        run = max(run, (m - r) * p[i])
        out[i] = min(1.0, run)
    return out


def merge(out):
    d = out / "iterations"
    files = [f for f in sorted(d.glob("iter*_detailed.csv")) if "_classical_" not in f.name]      # the classical-ML files of the earlier bake-off are not part of the paper
    if not files:
        sys.exit(f"no iteration files in {d}: run run_bakeoff.sh first")
    det = pd.concat([pd.read_csv(f) for f in files], ignore_index=True).sort_values(["iteration", "classifier"])
    pl = pd.concat([pd.read_csv(f) for f in sorted(d.glob("iter*_per_label.csv")) if "_classical_" not in f.name], ignore_index=True)
    counts = det.groupby("classifier").iteration.nunique()
    print(f"{det.iteration.nunique()} iterations; methods: {dict(counts)}")
    if det.iteration.nunique() < 16 or counts.nunique() > 1:
        print("WARNING: fewer than 16 iterations, or not every method has every iteration")
    cols = [m for m in METRICS if m in det.columns]
    s = det.groupby(["split_mode", "dataset", "classifier"])[cols].agg(["mean", "std"]).reset_index()
    s.columns = ["_".join(c).strip("_") for c in s.columns]
    det.to_csv(out / "sam_classifier_comparison_detailed.csv", index=False)
    pl.to_csv(out / "sam_classifier_comparison_per_label.csv", index=False)
    s.to_csv(out / "sam_classifier_comparison_summary.csv", index=False)
    pl.groupby(["classifier", "label"])[[m for m in METRICS if m in pl.columns]].mean().round(4).reset_index().to_csv(out / "sam_classifier_comparison_per_label_summary.csv", index=False)
    return det, pl


def paired_tests(det, out, n_boot=5000):
    w = det.pivot(index="iteration", columns="classifier", values="mcc")
    rng, rows = np.random.default_rng(0), []
    for m in [c for c in w.columns if c.startswith("MLP_")]:
        sub = []
        for b in [c for c in SEARCH if c in w.columns]:
            d = (w[m] - w[b]).values
            boots = [rng.choice(d, len(d)).mean() for _ in range(n_boot)]
            sub.append({"MLP": m, "baseline": b, "mlp_mcc": w[m].mean(), "baseline_mcc": w[b].mean(), "mean_diff": d.mean(), "ci95_lo": np.percentile(boots, 2.5),
                        "ci95_hi": np.percentile(boots, 97.5), "splits_mlp_wins": f"{int((d > 0).sum())}/{len(d)}", "p_wilcoxon": wilcoxon(d).pvalue if (d != 0).any() else 1.0})
        for r, hp in zip(sub, holm([r["p_wilcoxon"] for r in sub])):
            r["p_holm_3"] = hp
            rows.append(r)
    t = pd.DataFrame(rows)
    t.to_csv(out / "paired_tests_vs_baselines.csv", index=False)
    return t


def loss_curves(out):
    d = out / "iterations"
    cur = pd.concat([pd.read_csv(f) for f in sorted(d.glob("iter*_curves.csv.gz"))], ignore_index=True)
    tr = pd.concat([pd.read_csv(f) for f in sorted(d.glob("iter*_training.csv"))], ignore_index=True)
    tr.to_csv(out / "kept_epochs.csv", index=False)
    embs = [e for e in EMBS if e in set(cur.embedding)]

    def draw(ax, c, title, kept):
        va, trn = c.pivot(index="epoch", columns="iteration", values="val"), c.pivot(index="epoch", columns="iteration", values="train")
        for it in va.columns:
            ax.plot(va.index, va[it], color="#d62728", alpha=0.15, lw=0.6)
        ax.plot(trn.index, trn.mean(axis=1), color="#1f77b4", lw=1.6, label="train (mean of splits)")
        ax.plot(va.index, va.mean(axis=1), color="#d62728", lw=1.6, label="validation (mean of splits; thin = each split)")
        ax.plot(kept, np.full(len(kept), va.values.min() * 0.9), "|", color="k", ms=9, alpha=0.7)
        ax.axvline(np.median(kept), color="k", ls="--", lw=1.0)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(f"{title}: kept epoch median {np.median(kept):.0f} (range {kept.min()}-{kept.max()})", fontsize=9)

    fig, axs = plt.subplots(2, 3, figsize=(17, 8), squeeze=False)
    for ax, e in zip(axs.ravel(), embs):
        c = cur[cur.embedding == e].groupby(["iteration", "epoch"])[["train_loss", "val_loss"]].sum().reset_index().rename(columns={"train_loss": "train", "val_loss": "val"})
        draw(ax, c, e, tr[tr.embedding == e].best_epoch.values)
        ax.set_xlabel("epoch"); ax.set_ylabel("total loss (sum over labels)")
    axs[0, 0].legend(fontsize=7)
    fig.suptitle("BCE train / validation loss of the bake-off MLPs over all epochs; ticks = the epoch each split kept (early stopping), dashed = median", fontsize=11)
    fig.tight_layout(); fig.savefig(out / "train_val_loss_by_embedding.png", dpi=130); plt.close(fig)

    fig, axs = plt.subplots(len(embs), len(LABELS), figsize=(4.2 * len(LABELS), 2.6 * len(embs)), squeeze=False)
    for i, e in enumerate(embs):
        kept = tr[tr.embedding == e].best_epoch.values
        for j, lab in enumerate(LABELS):
            c = cur[(cur.embedding == e) & (cur.label == lab)].rename(columns={"train_loss": "train", "val_loss": "val"})
            draw(axs[i, j], c, f"{e} | {lab}", kept)
            axs[i, j].set_ylabel("loss") if j == 0 else None
            axs[i, j].set_xlabel("epoch") if i == len(embs) - 1 else None
    fig.suptitle("Per-label train / validation loss; ticks and dashed line = the epoch kept by early stopping on the TOTAL validation loss", fontsize=10)
    fig.tight_layout(); fig.savefig(out / "train_val_loss_by_label.png", dpi=110); plt.close(fig)
    return tr


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", help="the bake-off folder (default: $BAKEOFF_OUT, else <release>/benchmark)")
    a = ap.parse_args()
    import bakeoff
    out = Path(a.out or os.environ.get("BAKEOFF_OUT") or bakeoff.release_dir() / "benchmark")
    det, pl = merge(out)
    t = det.groupby("classifier")[[m for m in ("mcc", "mcc_valtuned") if m in det.columns]].agg(["mean", "std"]).round(3)
    print("\ntest MCC over the splits (mcc = the primary threshold: 0.5, HMM validation-tuned; mcc_valtuned = thresholds tuned on validation):")
    print(t.sort_values(("mcc", "mean"), ascending=False).to_string())
    print("\nper-label MCC, mean over the splits:")
    print(pl.pivot_table(index="classifier", columns="label", values="mcc")[LABELS].round(3).loc[t.sort_values(("mcc", "mean"), ascending=False).index].to_string())
    env = {**os.environ, "BAKEOFF_OUT": str(out), "BAKEOFF_CONTROL": bakeoff.CFG.get("control", "MLP_ProstT5")}
    for script in ("classifier_mcc_stats.py", "plot_cluster_mcc_zoom.py", "build_metrics_table.py", "benchmark_supplement_stats.py"):      # the last: the paired tests against the deployed model, the Friedman test of Fig 2A, the per-iteration numbers
        print(f"\n--- {script}")
        subprocess.run([sys.executable, script], cwd=HERE, env=env, check=True)
    p = paired_tests(det, out)
    print("\nMLP minus search baseline, paired over the splits (Holm over each MLP's three comparisons; the p-values are optimistic: the splits overlap):")
    fmt = p.assign(diff=p.mean_diff.map("{:+.3f}".format), ci=[f"[{lo:+.3f}, {hi:+.3f}]" for lo, hi in zip(p.ci95_lo, p.ci95_hi)], p=p.p_wilcoxon.map("{:.4f}".format), holm=p.p_holm_3.map("{:.4f}".format))
    print(fmt[["MLP", "baseline", "mlp_mcc", "baseline_mcc", "diff", "ci", "splits_mlp_wins", "p", "holm"]].round(3).to_string(index=False))
    tr = loss_curves(out)
    print("\nkept epoch (early stopping, patience %s), median [min-max] over the splits:" % bakeoff.PATIENCE)
    for e in [e for e in EMBS if e in set(tr.embedding)]:
        k = tr[tr.embedding == e]
        print(f"   {e:16s} {k.best_epoch.median():5.0f} [{k.best_epoch.min()}-{k.best_epoch.max()}]   early stopping fired in {k.stop_epoch.notna().mean():.0%} of the splits")
    print(f"\nwrote the merged csv files, the figures and paired_tests_vs_baselines.csv, kept_epochs.csv into {out}/")


if __name__ == "__main__":
    main()
