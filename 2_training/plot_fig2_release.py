#!/usr/bin/env python3
"""Figure 2 of a refresh release: A and B in one figure.

  python plot_fig2_release.py [--bakeoff <release>/benchmark] [--sam-vs-up <release>/sam_vs_up] [--out <png>]

A  cluster-split MCC of the six embedding MLPs and the three search-tool baselines (plot_cluster_mcc_zoom.draw_panel): Friedman box, stars = significantly worse than the deployed
   MLP (paired Wilcoxon, Holm), bracket = not significantly different from it.
B  MCC and AUPRC of the MLP trained on MTase proteins vs the one trained on the broad UP set, both scored on the held-out MTase test clusters (sam_vs_up_report.draw_arms): brackets =
   Wilcoxon over the splits, Holm-corrected within a panel.
No titles and no text under the panels (the description and the statistics are in the caption files written next to the panels); panel letters only.
Reads the CSVs the two reports wrote (bakeoff_report.py, sam_vs_up_report.py), so run it after both. Writes <release>/benchmark/Fig2.png and Fig2_caption.txt (the statistics sentences of both panels).
"""
import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bakeoff")
    ap.add_argument("--sam-vs-up")
    ap.add_argument("--out")
    a = ap.parse_args()
    import bakeoff
    rel = bakeoff.release_dir()
    bk = Path(a.bakeoff or os.environ.get("BAKEOFF_OUT") or rel / "benchmark")
    sv = Path(a.sam_vs_up or os.environ.get("SAM_VS_UP_OUT") or rel / "sam_vs_up")
    os.environ["BAKEOFF_OUT"] = str(bk)                                         # plot_cluster_mcc_zoom reads these when it is imported
    os.environ["BAKEOFF_CONTROL"] = bakeoff.CFG.get("control", "MLP_ProstT5")
    sys.path.insert(0, str(HERE))
    import plot_cluster_mcc_zoom as cz
    import sam_vs_up_report as su

    summ = pd.read_csv(sv / "sam_vs_up_summary.csv").set_index(["classifier", "label"])
    tests = pd.read_csv(sv / "sam_vs_up_paired_tests.csv")
    shown = [x for x in ("ProtT5_SAM", "ProtT5_UP") if (x, "Overall") in summ.index]

    fig = plt.figure(figsize=(13, 12.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 0.85], hspace=0.42, wspace=0.16, left=0.07, right=0.985, top=0.94, bottom=0.07)
    axA = fig.add_subplot(gs[0, :])
    cz.draw_panel(axA, "full", cz.REPRO_TOOLS)
    axA.legend(handles=cz.LEGEND, fontsize=10, loc="lower center", bbox_to_anchor=(0.5, 1.005), ncol=2, frameon=False)
    axB = [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    su.draw_arms(axB, summ, tests, shown)
    h, l = axB[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, fontsize=10, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.text(0.012, 0.975, "A", fontsize=26, fontweight="bold", va="top")
    fig.text(0.012, 0.425, "B", fontsize=26, fontweight="bold", va="top")
    cap = []      # the caption text for the draft: the two panels' sentences and statistics, written by the reports
    for f, tag in ((bk / "sam_classifier_cluster_mcc_zoom_caption.txt", "(A) "), (sv / "sam_vs_up_mcc_caption.txt", "(B) ")):
        if f.exists():
            cap.append(tag + f.read_text().strip())
    (bk / "Fig2_caption.txt").write_text("\n".join(cap) + "\n")
    out = Path(a.out) if a.out else bk / "Fig2.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
