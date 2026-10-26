#!/usr/bin/env python3
"""Table S2 -- Plant dark clusters (19), rendered as a figure (matches the
Table S3 visual convention: banded matplotlib table, no embedded title,
class-colored accent -- here the uniform "Other" substrate-class orange, since
every dark cluster is an Other-substrate, small-molecule family).

Two independent computational calls per cluster:
  - CLEAN's EC-centroid nearest-known-enzyme assignment (distance d;
    reliable/marginal/unreliable -> green/yellow/pink).
  - Horizyn's dominant reaction family (share of the cluster voting for that
    family, shaded green/yellow/pink by share; families flagged "NOT
    methylation" in bold red contain oxidoreductase/glycosyltransferase
    chemistry rather than SAM methylation).

Source: repo/6_dark_cluster_annotation/table2_dark_clusters.csv
Outputs (this dir): table_s2.csv, table_s2.png
"""
import textwrap
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = Path("/home/nroethler/Code/journal_article/repo/6_dark_cluster_annotation/table2_dark_clusters.csv")
HERE = Path(__file__).resolve().parent
LINEAGE = "Plant"
HEADER_COLOR = "#333333"
ACCENT = "#e6550d"   # "Other" substrate-class orange (Fig7 panel C dots)
GREEN, YELLOW, PINK = "#c7e9c0", "#fee391", "#fcbba1"

COLS = ["Cluster", "Nearest known enzyme\n(CLEAN)", "MTase?", "Confidence\n(CLEAN)",
        "Dominant reaction family (Horizyn)", "Family share, noise\n(Horizyn)", "Avg.\nHorizyn\nscore"]
WIDTHS = [0.075, 0.155, 0.045, 0.085, 0.37, 0.165, 0.055]
WRAP = {1: 26, 4: 58}


def tier(word):
    return {"reliable": GREEN, "marginal": YELLOW, "unreliable": PINK}.get(word, "white")


def share_tier(pct):
    if pct >= 70:
        return GREEN
    if pct >= 40:
        return YELLOW
    return PINK


def wrap(text, width):
    return "\n".join(textwrap.wrap(str(text), width=width, break_long_words=True,
                                    break_on_hyphens=True))


def main():
    df = pd.read_csv(SRC)
    df = df[df.lineage == LINEAGE].sort_values("n", ascending=False).reset_index(drop=True)

    rows = []
    for _, r in df.iterrows():
        conf_word = r["clean_confidence"].split()[0]
        share_pct = int(r["horizyn_share"].split("%")[0])
        not_meth = "NOT methylation" in r["horizyn_family"]
        rows.append({
            "Cluster": f"c{r['cluster']}\n(n={r['n']})",
            "Nearest known enzyme\n(CLEAN)": wrap(r["clean_enzyme"], WRAP[1]),
            "MTase?": r["clean_mtase"],
            "Confidence\n(CLEAN)": wrap(r["clean_confidence"], 14),
            "Dominant reaction family (Horizyn)": wrap(r["horizyn_family"], WRAP[4]),
            "Family share, noise\n(Horizyn)": wrap(r["horizyn_share"], 20),
            "Avg.\nHorizyn\nscore": f"{r['horizyn_avg_score']:.2f}",
            "_conf_color": tier(conf_word),
            "_share_color": share_tier(share_pct),
            "_not_meth": not_meth,
        })
    tbl = pd.DataFrame(rows)
    tbl.drop(columns=["_conf_color", "_share_color", "_not_meth"]).to_csv(
        HERE / "table_s2.csv", index=False)
    print(f"[wrote] table_s2.csv ({len(tbl)} rows)")

    # ---- rendered table ----
    nlines = lambda s: s.count("\n") + 1
    row_lines = [max(nlines(tbl.iloc[i][c]) for c in COLS) for i in range(len(tbl))]
    line_h = 0.34
    fig_h = sum(h * line_h + 0.12 for h in row_lines) + 0.6
    fig, ax = plt.subplots(figsize=(19, fig_h))
    ax.axis("off")
    t = ax.table(cellText=tbl[COLS].values, colLabels=COLS, cellLoc="left",
                 loc="center", bbox=[0, 0, 1, 1])
    t.auto_set_font_size(False); t.set_fontsize(9)
    for (r, c), cell in t.get_celld().items():
        cell.set_edgecolor("#dddddd")
        cell.set_width(WIDTHS[c])
        cell.PAD = 0.01
        if r == 0:
            cell.set_facecolor(HEADER_COLOR); cell.set_text_props(color="white", fontweight="bold")
            cell.set_height(0.08)
            continue
        row = tbl.iloc[r - 1]
        lines = row_lines[r - 1]
        cell.set_height(lines * line_h / fig_h + 0.02)
        col_name = COLS[c]
        if col_name == "Cluster":
            cell.set_facecolor(ACCENT); cell.set_text_props(color="white", fontweight="bold", va="center")
        elif col_name == "Confidence\n(CLEAN)":
            cell.set_facecolor(row["_conf_color"]); cell.set_text_props(va="center")
        elif col_name == "Family share, noise\n(Horizyn)":
            cell.set_facecolor(row["_share_color"]); cell.set_text_props(va="center")
        else:
            cell.set_facecolor("#ffffff" if (r % 2) else "#f4f4f4")
            cell.set_text_props(va="center")
            if col_name == "Dominant reaction family (Horizyn)" and row["_not_meth"]:
                cell.set_text_props(color="#b30000", fontweight="bold", va="center")
    fig.savefig(HERE / "table_s2.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("[wrote] table_s2.png")


if __name__ == "__main__":
    main()
