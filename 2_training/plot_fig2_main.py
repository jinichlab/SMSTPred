#!/usr/bin/env python3
"""Main Figure 2 — two panels in one figure:
   A  Representations & tools: cluster-disjoint MCC across the six embedding MLPs and
      three off-the-shelf search-tool baselines (from plot_cluster_mcc_zoom), with
      Friedman omnibus + Bonferroni-Dunn top-tier bracket.
   B  MTase-trained vs UP-trained ProstT5 MLP on held-out MTase proteins, MCC with Nadeau-Bengio
      corrected-paired-test significance (from plot_sam_vs_up_prostt5).

Everything else in the benchmark family is supplementary (MLP-vs-classical-ML zoom,
the MCC-vs-epoch and train/validation-loss grids, the MTase-vs-UP AUPRC panel, UP loss tuning).

Pure CSV -> PNG (no GPU). Reuses the panel draw functions from the two sibling
scripts so styling/statistics stay identical to their standalone versions.
"""
import os
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import plot_cluster_mcc_zoom as cz          # draw_panel(ax, ds); LEGEND
import plot_sam_vs_up_prostt5 as su         # draw_panel(ax, metric, summ, det, ylab)

HERE = Path(__file__).resolve().parent
# Optional re-run folders (environment variables); the defaults reproduce the committed Fig 2.
OUTA = HERE / os.environ.get("BAKEOFF_OUT", "sam_classifier_comparison_output")
OUTB = HERE / os.environ.get("SAM_VS_UP_OUT", "sam_vs_up_prostt5_comparison_output")


# Figure 2 style: four colours taken from recent Digital Discovery figures (the hue families that recur most across 455 figures in 88 articles are
# blue, orange, red and green; each colour below is the observed median shade of its family; see article/figures/Fig2_benchmark/palette_analysis/).
# Panel A uses blue + green and panel B orange + red, so the panels do not share a colour code and neither panel pairs red with green
# (indistinguishable for deuteranopes). Type is large because the figure is scaled to column width in print.
STYLE_A = dict(embed="#4578aa", search="#618d43", other="#969696", accent="#2e5a86",
               fs_tick=15, fs_ylabel=16, fs_value=13, fs_bracket=15, fs_stats=14, ytick_pad=4,
               value_bbox=True, stats_short=True, stats_xy=(0.99, 0.985))
STYLE_B = dict(colors={"SAM": "#ed8c44", "UP": "#c14445"}, fs_tick=16, fs_ylabel=16, fs_sig=15, fs_p=11.5, fs_legend=15,
               p_dy=0.085, sig_dy=0.008, ylim=1.40, bracket_cap=1.19)
FS_LEGEND = 15
FS_LETTER = 30


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None, help="output PNG (default: <bake-off output>/fig2_benchmark_main.png)")
    out = ap.parse_args().out
    summ = pd.read_csv(OUTB / "sam_vs_up_prostt5_comparison_summary.csv")
    det = pd.read_csv(OUTB / "sam_vs_up_prostt5_comparison_detailed.csv")

    fig = plt.figure(figsize=(12.5, 14.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.3, 0.85], hspace=0.62,
                          left=0.095, right=0.985, top=0.925, bottom=0.04)

    # ---- Panel A: representations & tools, cluster-disjoint MCC (full set) ----
    axA = fig.add_subplot(gs[0])
    cz.draw_panel(axA, "full", cz.REPRO_TOOLS, style=STYLE_A)
    axA.text(-0.085, 1.17, "A", transform=axA.transAxes, fontsize=FS_LETTER, fontweight="bold", va="top")
    axA.legend(handles=cz.legend_handles(STYLE_A), fontsize=FS_LEGEND, loc="lower center",
               bbox_to_anchor=(0.5, 1.005), ncol=2, frameon=True, handlelength=1.6, columnspacing=1.6)

    # ---- Panel B: MTase-trained vs UP-trained ProstT5 MLP, MCC ----
    axB = fig.add_subplot(gs[1])
    su.draw_panel(axB, "mcc", summ, det, "MCC (tuned threshold)", style=STYLE_B,
                  legend_kw=dict(loc="lower center", bbox_to_anchor=(0.5, 1.005), ncol=2, frameon=True, fontsize=FS_LEGEND,
                                 handlelength=1.6, columnspacing=1.6))
    axB.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    axB.text(-0.085, 1.27, "B", transform=axB.transAxes, fontsize=FS_LETTER, fontweight="bold", va="top")

    p = Path(out) if out else OUTA / "fig2_benchmark_main.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    print("wrote", p)


if __name__ == "__main__":
    main()
