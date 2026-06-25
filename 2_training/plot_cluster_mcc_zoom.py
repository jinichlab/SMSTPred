#!/usr/bin/env python3
"""Zoomed cluster-split MCC bar chart (full + balanced) with mean labels.

Reads sam_classifier_comparison_summary.csv and draws, per dataset, the
cluster-disjoint MCC for every classifier sorted high->low, with +/-1 SD error
bars, each bar annotated with its mean. Y-axis is zoomed to amplify the
differentiation between the top methods. Bars colored by family:
MLP (embedding models) / homology+profile baselines / classical ML.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

OUT = Path(__file__).resolve().parent / "sam_classifier_comparison_output"
s = pd.read_csv(OUT / "sam_classifier_comparison_summary.csv")

BASELINES = {"MMseqs2_1NN", "Foldseek_1NN", "Profile_HMM"}  # homology / profile


def color(clf):
    if clf.startswith("MLP"):
        return "#2c7fb8"          # embedding MLPs
    if clf in BASELINES:
        return "#e6550d"          # homology / profile baselines
    return "#969696"              # classical ML


YMIN = 0.50                       # zoom floor (lowest bar ~0.56)
datasets = ["full", "balanced"]
fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.4), sharey=True)

for ax, ds in zip(axes, datasets):
    sub = (s[(s.split_mode == "cluster") & (s.dataset == ds)]
           .sort_values("mcc_mean", ascending=False).reset_index(drop=True))
    x = range(len(sub))
    cols = [color(c) for c in sub.classifier]
    ax.bar(x, sub.mcc_mean, yerr=sub.mcc_std, capsize=3,
           color=cols, edgecolor="black", linewidth=0.6,
           error_kw=dict(ecolor="#444444", lw=1))
    for i, (mm, sdd) in enumerate(zip(sub.mcc_mean, sub.mcc_std)):
        ax.text(i, min(mm + sdd + 0.006, 0.998), f"{mm:.3f}",
                ha="center", va="bottom", fontsize=8, fontweight="bold")
    ax.set_xticks(list(x))
    ax.set_xticklabels(sub.classifier, rotation=45, ha="right", fontsize=8)
    ax.set_ylim(YMIN, 1.0)
    ax.set_title(f"cluster-disjoint / {ds}", fontsize=12, fontweight="bold")
    ax.set_ylabel("MCC (mean ± SD over 6 iters)")
    ax.grid(True, axis="y", alpha=0.3)
    ax.axhline(YMIN, color="black", lw=0.8)

legend = [mpatches.Patch(color="#2c7fb8", label="Embedding MLP"),
          mpatches.Patch(color="#e6550d", label="Homology / profile baseline"),
          mpatches.Patch(color="#969696", label="Classical ML")]
fig.legend(handles=legend, fontsize=8, loc="upper center", ncol=3,
           bbox_to_anchor=(0.5, 0.945), frameon=True)
fig.suptitle("SAM substrate-type classifiers — cluster-disjoint MCC (zoomed)",
             fontsize=13, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.90])
p = OUT / "sam_classifier_cluster_mcc_zoom.png"
fig.savefig(p, dpi=150)
print("wrote", p)
