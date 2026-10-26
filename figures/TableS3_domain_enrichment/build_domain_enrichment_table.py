#!/usr/bin/env python3
"""Table S3 — InterPro domain enrichment per predicted substrate class.

Source: repo/7_motif_analysis/domain_enrichment.csv (one-sided Fisher's exact
test, class-member InterPro domain presence/absence vs all other classes;
repo/7_motif_analysis/3_domain_enrichment.py). Reproduces the script's own
"most prevalent class-specific" cut: fold >= 2, top 12 by in-class protein
count per class, excluding the pan-methyltransferase SAM-binding superfamily
domain (IPR029063, fold = 1.0 by construction -- present in every class).

Outputs (this dir): domain_enrichment_table.csv, domain_enrichment_table.png
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = Path("/home/nroethler/Code/journal_article/repo/7_motif_analysis/domain_enrichment.csv")
HERE = Path(__file__).resolve().parent

CLASSES = ["DNA", "Protein", "RNA", "Other"]
COLORS = {"DNA": "#6a51a3", "RNA": "#cb181d", "Protein": "#238b45", "Other": "#e6550d"}
TOP_N = 12
PAN_MTASE_DOMAIN = "IPR029063"  # SAM-dependent MTase superfamily -- present in all classes


def fmt_p(p):
    return "<1e-300" if p == 0.0 else f"{p:.1e}"


def fmt_fold(f):
    return "inf" if np.isinf(f) else f"{f:.1f}"


def main():
    df = pd.read_csv(SRC)
    df = df[df["interpro"] != PAN_MTASE_DOMAIN]
    df = df[df["fold"] >= 2]
    df = df[df["name"].notna() & (df["name"].str.strip() != "")]

    out = []
    for cls in CLASSES:
        sub = df[df["class"] == cls].sort_values("n_in", ascending=False).head(TOP_N)
        for _, r in sub.iterrows():
            out.append({
                "Class": cls,
                "InterPro": r["interpro"],
                "Domain": r["name"],
                "Proteins": int(r["n_in"]),
                "% of class": f"{100 * r['pct_in']:.1f}",
                "Fold enrichment": fmt_fold(r["fold"]),
                "p-value": fmt_p(r["p"]),
            })
    tbl = pd.DataFrame(out)
    tbl.to_csv(HERE / "domain_enrichment_table.csv", index=False)
    print(tbl.to_string(index=False))
    print(f"\n[wrote] domain_enrichment_table.csv ({len(tbl)} rows)")

    # ---- rendered table ----
    cols = list(tbl.columns)
    fig, ax = plt.subplots(figsize=(17, 0.32 * len(tbl) + 0.3))
    ax.axis("off")
    t = ax.table(cellText=tbl.values, colLabels=cols, cellLoc="left",
                 loc="center", bbox=[0, 0, 1, 1])
    t.auto_set_font_size(False); t.set_fontsize(9.5)
    widths = {"Class": 0.06, "InterPro": 0.065, "Domain": 0.50, "Proteins": 0.075,
              "% of class": 0.085, "Fold enrichment": 0.10, "p-value": 0.105}
    for (r, c), cell in t.get_celld().items():
        cell.set_edgecolor("#dddddd")
        cell.set_width(widths[cols[c]])
        if r == 0:
            cell.set_facecolor("#333333"); cell.set_text_props(color="white", fontweight="bold")
        else:
            lab = tbl.iloc[r - 1]["Class"]
            if cols[c] == "Class":
                cell.set_facecolor(COLORS[lab]); cell.set_text_props(color="white", fontweight="bold")
            else:
                cell.set_facecolor("#ffffff" if (r % 2) else "#f4f4f4")
    fig.tight_layout()
    fig.savefig(HERE / "domain_enrichment_table.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("[wrote] domain_enrichment_table.png")


if __name__ == "__main__":
    main()
