#!/usr/bin/env python3
"""SAM-MTase count vs proteome size, Bacteria and Eukaryota separately, each
with a total-SAM panel and an Metabolite-substrate-SAM panel (4 panels total).
Actinomycetota highlighted in the Bacteria panels, Viridiplantae (plants) in
the Eukaryota panels, each against an OLS line of best fit over all points in
that domain -- to see whether the highlighted lineage sits on, above, or below
the overall scaling trend.

Input:  proteome_vs_sam_bylabel.csv (organism, domain, proteome_size, n_sam,
        n_DNA, n_RNA, n_Protein, n_Metabolite -- from the existing per-organism
        tree-of-life sweep), selected_1259_organisms.csv (Taxonomic lineage,
        joined by organism-name stem)
Output: proteome_vs_sam_domain_highlight.png
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IO = Path(__file__).resolve().parent
META_CSV = Path(os.environ.get("SMSTPRED_DATA", Path(__file__).resolve().parents[2] / "data")) / "selected_1259_organisms.csv"

BASE_COLOR = "#707070"
ACTINO_COLOR = "#ed8c44"
PLANT_COLOR = "#618d43"
FIT_COLOR = "#222222"


def load():
    d = pd.read_csv(IO / "proteome_vs_sam_bylabel.csv")
    meta = pd.read_csv(META_CSV)
    meta["stem"] = meta["Organism"].astype(str).str.replace(" ", "_", regex=False).str.replace("/", "_", regex=False)
    d = d.merge(meta[["stem", "Taxonomic lineage"]], left_on="organism", right_on="stem", how="left")
    lin = d["Taxonomic lineage"].fillna("")
    d["is_actino"] = lin.str.contains("Actinomycetota")
    d["is_plant"] = lin.str.contains("Viridiplantae")
    return d


def panel(ax, sub, ycol, highlight_mask, highlight_color, highlight_label, other_label, ylabel):
    x = sub["proteome_size"].to_numpy(float)
    y = sub[ycol].to_numpy(float)
    hi = highlight_mask.to_numpy()

    ax.scatter(x[~hi], y[~hi], s=16, c=BASE_COLOR, alpha=0.5, linewidths=0,
               zorder=2, label=f"{other_label} (n={int((~hi).sum())})")
    ax.scatter(x[hi], y[hi], s=34, c=highlight_color, alpha=0.85,
               edgecolors="black", linewidths=0.4, marker="^", zorder=4,
               label=f"{highlight_label} (n={int(hi.sum())})")

    # separate OLS fit per group, so the two groups' scaling rates can be
    # compared directly rather than only against a single pooled trend
    slope_hi, ic_hi = np.polyfit(x[hi], y[hi], 1)
    slope_rest, ic_rest = np.polyfit(x[~hi], y[~hi], 1)
    xs_hi = np.array([x[hi].min(), x[hi].max()])
    xs_rest = np.array([x[~hi].min(), x[~hi].max()])
    ax.plot(xs_hi, slope_hi * xs_hi + ic_hi, "-", color=highlight_color, lw=2.0,
            zorder=3, label=f"{highlight_label} fit")
    ax.plot(xs_rest, slope_rest * xs_rest + ic_rest, "--", color=FIT_COLOR, lw=1.6,
            zorder=3, label=f"{other_label} fit")

    # slope as a memorable rate -- SAM-MTases per 100 proteins -- compared
    # between the highlighted lineage and the rest of the domain
    per100_hi = slope_hi * 100
    per100_rest = slope_rest * 100
    ratio = per100_hi / per100_rest if per100_rest else float("nan")
    box = dict(boxstyle="round,pad=0.35", fc="white", ec="#cccccc", alpha=0.9)
    ax.text(0.035, 0.965,
            f"{highlight_label}: {per100_hi:.2f} / 100 proteins\n"
            f"{other_label}: {per100_rest:.2f} / 100 proteins\n"
            f"({ratio:.1f}× higher)",
            transform=ax.transAxes, va="top", ha="left", fontsize=10,
            color=highlight_color, fontweight="bold", bbox=box)

    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_xlim(0, x.max() * 1.05)
    ax.set_ylim(0, y.max() * 1.1)
    ax.grid(True, ls=":", alpha=0.35)
    ax.tick_params(labelsize=9)
    ax.legend(loc="lower right", fontsize=10, frameon=True)


def main():
    d = load()
    bact = d[d.domain == "Bacteria"].copy()
    euk = d[d.domain == "Eukaryota"].copy()

    fig, axes = plt.subplots(2, 2, figsize=(13, 10.5))

    panel(axes[0, 0], bact, "n_sam", bact["is_actino"], ACTINO_COLOR, "Actinomycetota",
          "Other bacteria", "Total SAM MTases")

    panel(axes[0, 1], bact, "n_Metabolite", bact["is_actino"], ACTINO_COLOR, "Actinomycetota",
          "Other bacteria", "Metabolite-substrate SAM MTases")

    panel(axes[1, 0], euk, "n_sam", euk["is_plant"], PLANT_COLOR, "Viridiplantae",
          "Other eukaryotes", "Total SAM MTases")

    panel(axes[1, 1], euk, "n_Metabolite", euk["is_plant"], PLANT_COLOR, "Viridiplantae",
          "Other eukaryotes", "Metabolite-substrate SAM MTases")

    for letter, ax in zip("ABCD", axes.flat):
        ax.text(0.0, 1.02, letter, transform=ax.transAxes, fontsize=20,
                fontweight="bold", ha="left", va="bottom")

    for ax in axes[1, :]:
        ax.set_xlabel("Proteome size (number of proteins)", fontsize=11)

    fig.tight_layout()
    fig.savefig(IO / "proteome_vs_sam_domain_highlight.png", dpi=200, bbox_inches="tight")
    print("wrote proteome_vs_sam_domain_highlight.png")


if __name__ == "__main__":
    main()
