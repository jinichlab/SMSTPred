#!/usr/bin/env python3
"""Main Figure 2 — two panels in one figure:
   A  Representations & tools: cluster-disjoint MCC across the six embedding MLPs and
      three off-the-shelf search-tool baselines (from plot_cluster_mcc_zoom), with
      Friedman omnibus + Bonferroni-Dunn top-tier bracket.
   B  SAM-trained vs UP-trained ProstT5 MLP on held-out SAM, MCC with Nadeau-Bengio
      corrected-paired-test significance (from plot_sam_vs_up_prostt5).

Everything else in the benchmark family is supplementary (MLP-vs-classical-ML zoom,
balanced classifier zoom, full/balanced MCC-vs-epoch grids, the SAM-vs-UP AUPRC
panel, UP loss tuning).

Pure CSV -> PNG (no GPU). Reuses the panel draw functions from the two sibling
scripts so styling/statistics stay identical to their standalone versions.
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import plot_cluster_mcc_zoom as cz          # draw_panel(ax, ds); LEGEND
import plot_sam_vs_up_prostt5 as su         # draw_panel(ax, metric, summ, det, ylab)

HERE = Path(__file__).resolve().parent
OUTA = HERE / "sam_classifier_comparison_output"
OUTB = HERE / "sam_vs_up_prostt5_comparison_output"


def main():
    summ = pd.read_csv(OUTB / "sam_vs_up_prostt5_comparison_summary.csv")
    det = pd.read_csv(OUTB / "sam_vs_up_prostt5_comparison_detailed.csv")

    fig = plt.figure(figsize=(13.5, 11.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.18, 0.86], hspace=0.62,
                          left=0.065, right=0.985, top=0.90, bottom=0.05)

    # ---- Panel A: representations & tools, cluster-disjoint MCC (full set) ----
    axA = fig.add_subplot(gs[0])
    cz.draw_panel(axA, "full", cz.REPRO_TOOLS)
    axA.text(-0.055, 1.20, "A", transform=axA.transAxes, fontsize=19, fontweight="bold", va="top")
    axA.legend(handles=cz.LEGEND, fontsize=9, loc="lower center",
               bbox_to_anchor=(0.5, 1.005), ncol=len(cz.LEGEND), frameon=True)

    # ---- Panel B: SAM-trained vs UP-trained ProstT5 MLP, MCC ----
    axB = fig.add_subplot(gs[1])
    su.draw_panel(axB, "mcc", summ, det, "MCC (tuned threshold)")
    axB.text(-0.055, 1.14, "B", transform=axB.transAxes, fontsize=19, fontweight="bold", va="top")

    p = OUTA / "fig2_benchmark_main.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    print("wrote", p)


if __name__ == "__main__":
    main()
