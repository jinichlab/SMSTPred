#!/usr/bin/env python3
"""
Stage 2 -- the statistics and figures of the SAM-vs-UP comparison (after run_sam_vs_up.sh):

  python sam_vs_up_report.py [--out <release>/sam_vs_up]

  sam_vs_up_{detailed,per_label,summary}.csv   the 16 splits merged; mean +- SD of every arm and label
  sam_vs_up_paired_tests.csv                   per label and overall: UP minus SAM, SAMclean minus SAM, UP minus SAMclean -- mean difference, splits won, bootstrap interval, three tests:
                                                 * Wilcoxon signed-rank over the splits,
                                                 * the paired t-test (assumes independent splits: too optimistic here),
                                                 * the Nadeau-Bengio corrected resampled t-test: the variance is inflated by (1/n + n_test/n_train) because the training sets of repeated
                                                   splits overlap. n_test/n_train is taken from the SAM arm (the smaller training set, so the conservative choice); the splits are drawn
                                                   as whole clusters, which the formula does not model, so read it as an approximation, not as exact.
  sam_vs_up_mcc.png                            MCC (and AUPRC) per label for the arms, mean +- SD over the splits
  up_train_val_loss.png, kept_epochs.csv       train / validation loss of the UP and SAMclean MLPs over all epochs and the epoch early stopping kept
All tests are two-sided. The arms and the design: sam_vs_up.py.
"""
import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

LABELS = ["Metabolite", "Protein", "RNA", "DNA"]
ARMS = ["ProtT5_SAM", "ProtT5_UP", "ProtT5_SAMclean"]
COLORS = {"ProtT5_SAM": "#4578aa", "ProtT5_UP": "#618d43", "ProtT5_SAMclean": "#7aa6c8"}


def holm_adjust(p):
    """Holm step-down adjusted p-values (same order as the input)."""
    p = np.asarray(p, float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (len(p) - rank) * p[i]))
        adj[i] = run
    return adj


def draw_arms(axs, summ, t, shown):
    """The two panels (MCC, AUPRC) of the MTase-trained vs UP-trained comparison into axs; summ = sam_vs_up_summary.csv indexed by (classifier, label), t = sam_vs_up_paired_tests.csv.
    Brackets: Wilcoxon signed-rank over the splits, Holm-corrected across the five comparisons of a panel. Returns the Nadeau-Bengio p-values of the comparisons shown."""
    keys = ["Overall"] + LABELS
    nb_p = []
    for ax, metric, title in zip(axs, ("mcc", "auprc"), ("MCC (fixed 0.5 threshold)", "AUPRC")):
        w = 0.8 / len(shown)
        top = []
        for i, arm in enumerate(shown):
            m = [summ.loc[(arm, lab), f"{metric}_mean"] for lab in keys]
            s = [summ.loc[(arm, lab), f"{metric}_std"] for lab in keys]
            top.append(np.array(m) + np.array(s))
            ax.bar(np.arange(5) + (i - (len(shown) - 1) / 2) * w, m, w, yerr=s, capsize=2, color=COLORS[arm], edgecolor="k", linewidth=0.5,
                   label={"ProtT5_SAM": "MTase-trained", "ProtT5_UP": "UP-trained (homologues of the MTase test clusters removed)"}[arm])
        if len(shown) == 2:       # Wilcoxon signed-rank over the splits, Holm-corrected across the five comparisons of this panel
            tt = t[(t.comparison == "ProtT5_UP minus ProtT5_SAM") & (t.metric == metric)].set_index("label").loc[keys]
            padj = holm_adjust(tt.p_wilcoxon.values)
            nb_p += list(tt.p_nadeau_bengio.values)
            for k, (lab, pa) in enumerate(zip(keys, padj)):
                y = min(max(top[0][k], top[1][k]) + 0.012, 1.045)
                ax.plot([k - w / 2, k - w / 2, k + w / 2, k + w / 2], [y, y + 0.008, y + 0.008, y], color="k", lw=0.9)
                ax.text(k, y + 0.010, "ns" if pa >= 0.05 else "**" if pa < 0.01 else "*", ha="center", va="bottom", fontsize=9)
                if pa < 0.05:
                    print(f"   {metric} {lab}: UP minus MTase-trained {tt.loc[lab, 'mean_diff']:+.3f}, Wilcoxon p = {tt.loc[lab, 'p_wilcoxon']:.4f}, Holm-adjusted {pa:.3f}")
        ax.set_xticks(np.arange(5)); ax.set_xticklabels(keys); ax.set_ylim(0.5, 1.09); ax.set_yticks(np.arange(0.5, 1.001, 0.1)); ax.set_title(title); ax.set_ylabel(title.split(" (")[0])
    return nb_p


def nadeau_bengio(d, rho):
    """Two-sided p of the corrected resampled t-test on the paired differences d (n splits); rho = n_test / n_train."""
    n = len(d)
    sd = d.std(ddof=1)
    if sd == 0:
        return np.nan, 1.0
    t = d.mean() / np.sqrt((1 / n + rho) * sd ** 2)
    return t, 2 * stats.t.sf(abs(t), n - 1)


def paired_row(a, b, rho, rng, label, metric, na, nb):
    d = (a - b).values
    boots = [rng.choice(d, len(d)).mean() for _ in range(5000)]
    t_naive, p_naive = stats.ttest_rel(a.values, b.values)
    t_nb, p_nb = nadeau_bengio(d, rho)
    return {"comparison": f"{na} minus {nb}", "label": label, "metric": metric, "mean_a": a.mean(), "mean_b": b.mean(), "mean_diff": d.mean(), "ci95_lo": np.percentile(boots, 2.5),
            "ci95_hi": np.percentile(boots, 97.5), "splits_a_wins": f"{int((d > 0).sum())}/{len(d)}", "p_wilcoxon": stats.wilcoxon(d).pvalue if (d != 0).any() else 1.0,
            "p_paired_t_naive": p_naive, "p_nadeau_bengio": p_nb}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out")
    a = ap.parse_args()
    import bakeoff
    rel = bakeoff.release_dir()
    out = Path(a.out or os.environ.get("SAM_VS_UP_OUT") or rel / "sam_vs_up")
    d = out / "iterations"
    det = pd.concat([pd.read_csv(f) for f in sorted(d.glob("iter*_detailed.csv"))], ignore_index=True)
    pl = pd.concat([pd.read_csv(f) for f in sorted(d.glob("iter*_per_label.csv"))], ignore_index=True)
    arms = [x for x in ARMS if x in set(det.classifier)]
    shown = [x for x in ("ProtT5_SAM", "ProtT5_UP") if x in arms]       # the figure shows the MTase-trained and the UP-trained arm only (the homologue-removal sensitivity arm is in the tables)
    n_it = det.iteration.nunique()
    print(f"{n_it} splits; arms {arms}")
    det.to_csv(out / "sam_vs_up_detailed.csv", index=False)
    pl.to_csv(out / "sam_vs_up_per_label.csv", index=False)
    summ = pd.concat([det.assign(label="Overall").drop(columns=["split_mode", "dataset"]), pl.drop(columns=["split_mode", "dataset"])]).groupby(["classifier", "label"])[["mcc", "auprc", "auroc", "f1", "accuracy"]].agg(["mean", "std"])
    summ.columns = ["_".join(c) for c in summ.columns]
    summ.reset_index().to_csv(out / "sam_vs_up_summary.csv", index=False)
    print("\nmean +- SD over the splits (MCC at a fixed 0.5 threshold; Overall = the macro mean over the four labels):")
    print(pd.DataFrame({arm: {lab: f"{summ.loc[(arm, lab), 'mcc_mean']:.3f} +- {summ.loc[(arm, lab), 'mcc_std']:.3f}" for lab in ["Overall"] + LABELS} for arm in arms}).to_string())

    sam_det = det[det.classifier == "ProtT5_SAM"] if "ProtT5_SAM" in arms else det
    rho = float(sam_det.n_test.mean() / sam_det.n_train.mean())
    print(f"\nNadeau-Bengio rho = n_test / n_train = {rho:.3f} (SAM arm: mean test {sam_det.n_test.mean():.0f}, mean train {sam_det.n_train.mean():.0f})")
    rng = np.random.default_rng(0)
    rows = []
    pairs = [(x, y) for x, y in (("ProtT5_UP", "ProtT5_SAM"), ("ProtT5_SAMclean", "ProtT5_SAM"), ("ProtT5_UP", "ProtT5_SAMclean")) if x in arms and y in arms]
    for na, nb in pairs:
        for lab in ["Overall"] + LABELS:
            for metric in ("mcc", "auprc"):
                if lab == "Overall":
                    pa, pb = det[det.classifier == na].set_index("iteration")[metric], det[det.classifier == nb].set_index("iteration")[metric]
                else:
                    pa = pl[(pl.classifier == na) & (pl.label == lab)].set_index("iteration")[metric]
                    pb = pl[(pl.classifier == nb) & (pl.label == lab)].set_index("iteration")[metric]
                pa, pb = pa.dropna(), pb.dropna()
                common = pa.index.intersection(pb.index)
                rows.append(paired_row(pa.loc[common], pb.loc[common], rho, rng, lab, metric, na, nb))
    t = pd.DataFrame(rows)
    t.to_csv(out / "sam_vs_up_paired_tests.csv", index=False)
    fp = lambda x: "<0.001" if x < 0.001 else f"{x:.3f}"
    for na, nb in pairs:
        s = t[(t.comparison == f"{na} minus {nb}") & (t.metric == "mcc")]
        print(f"\n{na} minus {nb}, MCC, paired over the splits (positive = the first arm is better):")
        print(pd.DataFrame({"label": s.label, "A": s.mean_a.round(3), "B": s.mean_b.round(3), "diff": s.mean_diff.map("{:+.3f}".format),
                            "95% CI": [f"[{lo:+.3f}, {hi:+.3f}]" for lo, hi in zip(s.ci95_lo, s.ci95_hi)], "A wins": s.splits_a_wins, "Wilcoxon": s.p_wilcoxon.map(fp),
                            "paired t": s.p_paired_t_naive.map(fp), "Nadeau-Bengio": s.p_nadeau_bengio.map(fp)}).to_string(index=False))

    # ---- figure: per-label MCC and AUPRC of the MTase-trained and the UP-trained arm, with the paired tests of UP minus MTase
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.0))
    nb_p = draw_arms(axs, summ, t, shown)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, fontsize=9, frameon=False)       # under the panels, clear of the bars
    note = ("Paired over the splits: Wilcoxon signed-rank test, Holm-corrected across the five comparisons of a panel (* p < 0.05, ** p < 0.01, ns = not significant)."
            + (f"\nNadeau-Bengio corrected t-test (the splits overlap): p ≥ {min(nb_p):.2f} for every comparison shown." if nb_p else ""))
    # no title and no text under the figure (they belong in the caption): the description and the statistics go to a caption file
    (out / "sam_vs_up_mcc_caption.txt").write_text(f"ProtT5 MLP trained on MTase proteins vs on UP (Swiss-Prot enzymes), scored on the same held-out MTase test clusters (mean ± SD over {n_it} cluster-disjoint splits). " + note.replace("\n", " ") + "\n")
    fig.tight_layout(rect=[0, 0.06, 1, 1]); fig.savefig(out / "sam_vs_up_mcc.png", dpi=130); plt.close(fig)

    # ---- UP / SAMclean training curves and kept epochs
    cf = sorted(d.glob("iter*_curves.csv.gz"))
    if cf:
        cur = pd.concat([pd.read_csv(f) for f in cf], ignore_index=True)
        tr = pd.concat([pd.read_csv(f) for f in sorted(d.glob("iter*_training.csv"))], ignore_index=True)
        tr.to_csv(out / "kept_epochs.csv", index=False)
        trained = [x for x in ("ProtT5_UP", "ProtT5_SAMclean") if x in set(cur.arm)]
        fig, axs = plt.subplots(1, len(trained), figsize=(7.5 * len(trained), 4.5), squeeze=False)
        for ax, arm in zip(axs[0], trained):
            c = cur[cur.arm == arm].groupby(["iteration", "epoch"])[["train_loss", "val_loss"]].sum().reset_index()
            va, trn = c.pivot(index="epoch", columns="iteration", values="val_loss"), c.pivot(index="epoch", columns="iteration", values="train_loss")
            for it in va.columns:
                ax.plot(va.index, va[it], color="#d62728", alpha=0.15, lw=0.6)
            ax.plot(trn.index, trn.mean(axis=1), color="#1f77b4", lw=1.6, label="train (mean of splits)")
            ax.plot(va.index, va.mean(axis=1), color="#d62728", lw=1.6, label="validation (mean of splits; thin = each split)")
            kept = tr[tr.arm == arm].best_epoch.values
            ax.plot(kept, np.full(len(kept), va.values.min() * 0.9), "|", color="k", ms=9, alpha=0.7); ax.axvline(np.median(kept), color="k", ls="--", lw=1)
            ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("epoch"); ax.set_ylabel("total loss (sum over labels)")
            ax.set_title(f"{arm}: kept epoch median {np.median(kept):.0f} (range {kept.min()}-{kept.max()})", fontsize=9)
        axs[0, 0].legend(fontsize=7)
        fig.tight_layout(); fig.savefig(out / "up_train_val_loss.png", dpi=130); plt.close(fig)
        print("\nkept epoch (early stopping), median [min-max] over the splits:")
        for arm in trained:
            k = tr[tr.arm == arm]
            print(f"   {arm:16s} {k.best_epoch.median():5.0f} [{k.best_epoch.min()}-{k.best_epoch.max()}]  train proteins {int(k.n_train.mean())}, validation {int(k.n_val.mean())}; early stopping fired in {k.stop_epoch.notna().mean():.0%}")
    print(f"\nwrote the merged csv files, sam_vs_up_paired_tests.csv, the figures into {out}/")


if __name__ == "__main__":
    main()
