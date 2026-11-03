#!/usr/bin/env python3
"""Table 1 -- cluster-disjoint benchmark, all metrics (MCC/AUPRC/AUROC/F1/accuracy)
for the nine representation and search-tool methods (Fig 2A's method set),
full training set, mean +/- SD over sixteen iterations. No embedded title (the
manuscript caption carries the description); family-shaded rows (embedding MLPs blue,
search tools green), the best value of each metric highlighted, no bold row.

Source: sam_classifier_comparison_output/sam_classifier_comparison_summary.csv
        (written by sam_classifier_comparison.py)
Output: sam_classifier_comparison_output/metrics_table.png
"""
import os
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE / os.environ.get("BAKEOFF_OUT", "sam_classifier_comparison_output")   # BAKEOFF_OUT: optional re-run folder
SRC = OUT_DIR / "sam_classifier_comparison_summary.csv"

CONTROL = "MLP_ProstT5"
EMBED_MLPS = {"MLP_ProstT5", "MLP_SaProt", "MLP_ProtT5", "MLP_ESM-IF",
              "MLP_ProtT5+ESM-IF", "MLP_ProstT5+ESM-IF"}
SEARCH_TOOLS = {"MMseqs2_1NN", "Profile_HMM", "Foldseek_1NN"}
METHODS = EMBED_MLPS | SEARCH_TOOLS
METRICS = ["mcc", "auprc", "auroc", "f1", "accuracy"]
METRIC_LABELS = {"mcc": "MCC", "auprc": "AUPRC", "auroc": "AUROC", "f1": "F1", "accuracy": "Accuracy"}


def main():
    df = pd.read_csv(SRC)
    sub = df[(df.split_mode == "cluster") & (df.dataset == "full") & (df.classifier.isin(METHODS))].copy()
    sub = sub.sort_values("mcc_mean", ascending=False).reset_index(drop=True)

    cols = ["Classifier"] + [METRIC_LABELS[m] for m in METRICS]
    rows_fmt = []
    best_col = {m: sub[f"{m}_mean"].idxmax() for m in METRICS}
    for _, r in sub.iterrows():
        row = [r["classifier"]] + [f"{r[f'{m}_mean']:.3f} ± {r[f'{m}_std']:.3f}" for m in METRICS]
        rows_fmt.append(row)

    fig, ax = plt.subplots(figsize=(13.5, 0.55 * len(sub) + 0.3))
    ax.axis("off")
    t = ax.table(cellText=rows_fmt, colLabels=cols, cellLoc="center", loc="center", bbox=[0, 0, 1, 1])
    t.auto_set_font_size(False); t.set_fontsize(11)
    widths = {"Classifier": 0.22}
    other_w = (1 - widths["Classifier"]) / len(METRICS)
    for m in METRICS:
        widths[METRIC_LABELS[m]] = other_w

    for (r, c), cell in t.get_celld().items():
        cell.set_edgecolor("#dddddd")
        cell.set_width(widths[cols[c]])
        if r == 0:
            cell.set_facecolor("#333333"); cell.set_text_props(color="white", fontweight="bold")
            continue
        clf = sub.iloc[r - 1]["classifier"]
        is_best = c >= 1 and (r - 1) == best_col[METRICS[c - 1]]
        cell.set_facecolor("#fff2a8" if is_best else "#d3e0ee" if clf in EMBED_MLPS else "#dce8d3")
        if c == 0:
            cell.set_text_props(ha="left")
            cell.PAD = 0.02

    fig.tight_layout()
    out = OUT_DIR / "metrics_table.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
