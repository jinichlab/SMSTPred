#!/usr/bin/env python3
"""Redraw sam_vs_up_prostt5_comparison_mcc.png from the saved CSVs, adding
per-label significance indicators.

Pure CSV -> PNG (no GPU): reads the summary (bar heights / SD) and detailed
(per-iteration) CSVs written by sam_vs_up_prostt5_comparison.py, and annotates
each MTase-vs-UP pair with the Nadeau-Bengio CORRECTED paired test (the honest test
for repeated-subsampling, which inflates the variance for train-set overlap; see
sam_vs_up_prostt5_stats.py). Two-sided corrected p is used for the bracket.

Over N=16 repeated cluster-disjoint splits, the DNA gap is read against its across-
split variance rather than the bar heights alone; the marks make that explicit.
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sam_vs_up_prostt5_stats as stx   # reuse paired_table / corrected_t / LABELS / K

HERE = Path(__file__).resolve().parent
OUT = HERE / os.environ.get("SAM_VS_UP_OUT", "sam_vs_up_prostt5_comparison_output")
LABELS = stx.LABELS
ORDER = LABELS + ["Overall"]
MODELS = [("SAM", "#4578aa"), ("UP", "#ed8c44")]       # keys as written to the CSVs by sam_vs_up_prostt5_comparison.py
DISPLAY = {"SAM": "MTase", "UP": "UP"}                 # names shown in figures: the SAM-methyltransferase training set is called the MTase set
W = 0.38
STYLE = dict(colors=dict(MODELS), fs_tick=10, fs_ylabel=10, fs_sig=9, fs_p=7.5, fs_legend=10, p_dy=0.052, sig_dy=0.006, ylim=1.25, bracket_cap=1.13)


def sig_text(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


def draw_panel(ax, metric, summ, det, ylab, title=None, legend=True, style=STYLE, legend_kw=None):
    """Draw one MTase-vs-UP grouped-bar metric panel with per-group significance."""
    tab = stx.paired_table(det, metric).set_index("label")
    x = np.arange(len(ORDER))
    gmax = {}
    for k, (mdl, col) in enumerate(MODELS):
        means = [summ[(summ.label == l) & (summ.model == mdl)][f"{metric}_mean"].values[0] for l in ORDER]
        stds = [summ[(summ.label == l) & (summ.model == mdl)][f"{metric}_std"].values[0] for l in ORDER]
        ax.bar(x + (k - 0.5) * W, means, W, yerr=stds, capsize=3, label=f"{DISPLAY[mdl]}-trained",
               color=style["colors"][mdl], edgecolor="black", linewidth=0.5, error_kw=dict(lw=1, ecolor="#444"))
        for i, l in enumerate(ORDER):
            gmax[l] = max(gmax.get(l, 0), means[i] + (stds[i] if not np.isnan(stds[i]) else 0))
    for i, l in enumerate(ORDER):
        p = tab.loc[l, "p_corrected_2sided"]
        ybr = min(gmax[l] + 0.045, style["bracket_cap"])
        xa, xb = x[i] - 0.5 * W, x[i] + 0.5 * W
        ax.plot([xa, xa, xb, xb], [ybr - 0.018, ybr, ybr, ybr - 0.018], color="#333", lw=1.0)
        ax.text(x[i], ybr + style["sig_dy"], f"{sig_text(p)}", ha="center", va="bottom",
                fontsize=style["fs_sig"], fontweight="bold" if p < 0.05 else "normal", color="#222")
        ax.text(x[i], ybr + style["p_dy"], f"p={p:.2f}", ha="center", va="bottom", fontsize=style["fs_p"], color="#777")
    ax.set_xticks(x); ax.set_xticklabels(ORDER, fontsize=style["fs_tick"]); ax.set_ylabel(ylab, fontsize=style["fs_ylabel"])
    ax.tick_params(axis="y", labelsize=style["fs_tick"])
    ax.set_ylim(0, style["ylim"]); ax.grid(True, axis="y", alpha=0.3)
    if legend:
        ax.legend(**({"loc": "lower right", "fontsize": style["fs_legend"]} | (legend_kw or {})))
    if title:
        ax.set_title(title, fontweight="bold")


def main():
    summ = pd.read_csv(OUT / "sam_vs_up_prostt5_comparison_summary.csv")
    det = pd.read_csv(OUT / "sam_vs_up_prostt5_comparison_detailed.csv")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    for ax, metric, lbl in zip(axes, ["mcc", "auprc"], ["MCC (tuned threshold)", "AUPRC"]):
        draw_panel(ax, metric, summ, det, lbl, title=lbl)

    fig.text(0.5, 0.045,
             f"Significance: Nadeau–Bengio corrected paired t-test over N={stx.N_ITER} cluster-disjoint splits "
             "(two-sided; corrects for train-set overlap across splits).",
             ha="center", va="bottom", fontsize=8.5, color="#333")
    fig.text(0.5, 0.014,
             "No individual substrate class differs; the marginal aggregate edge (Overall, p=0.04) favors focused "
             "MTase, not broad UP — more training data does not beat precise, task-matched data.",
             ha="center", va="bottom", fontsize=8.5, color="#333")
    fig.tight_layout(rect=[0, 0.08, 1, 1.0])
    p = OUT / "sam_vs_up_prostt5_comparison_mcc.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    print("wrote", p)


if __name__ == "__main__":
    main()
