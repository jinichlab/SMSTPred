#!/usr/bin/env python3
"""Static global UMAP panel (Fig 4) from the precomputed coords CSV.

Lightweight re-render (base env, no datamapplot): reads the final per-point
Substrate call straight from sam_mtase_umap_prostt5_coords.csv and draws the same
substrate-colored scatter as substrate_umap_datamap.create_static_plot, but with a
title that states the takeaway (distinct, substrate-coherent islands) and a pointer
to the interactive version in the repo. Marker sizes / alpha / draw order unchanged."""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
COORDS = HERE / "sam_mtase_umap_prostt5_coords.csv"
OUT = HERE / "sam_mtase_umap_prostt5_static_datamap.png"

SUBSTRATE_COLORS = {"RNA": "#cb181d", "DNA": "#6a51a3", "Protein": "#238b45",
                    "Other": "#e6550d", "Mixed": "#636363", "None": "#cccccc"}
DRAW_ORDER = ["None", "Mixed", "Other", "Protein", "DNA", "RNA"]
DIM = {"None": (0.25, 5), "Mixed": (0.40, 6)}     # (alpha, size) for de-emphasized classes


def main():
    df = pd.read_csv(COORDS)
    df["Substrate"] = df["Substrate"].replace({"Unassigned": "None"})

    fig, ax = plt.subplots(figsize=(14, 10))
    handles = {}
    for k, sub in enumerate(DRAW_ORDER):
        m = df["Substrate"] == sub
        if m.sum() == 0:
            continue
        alpha, size = DIM.get(sub, (0.7, 8))
        handles[sub] = ax.scatter(df.loc[m, "UMAP_1"], df.loc[m, "UMAP_2"],
                                  c=SUBSTRATE_COLORS[sub], label=f"{sub} (n={int(m.sum()):,})",
                                  alpha=alpha, s=size, edgecolors="none", zorder=2 + k)
    ax.set_xlabel("UMAP 1", fontsize=12); ax.set_ylabel("UMAP 2", fontsize=12)
    ax.set_title("Tree-of-life SAM-MTases resolve into distinct, substrate-coherent islands",
                 fontsize=15, fontweight="bold", pad=26)
    ax.text(0.5, 1.012,
            f"UMAP of {len(df):,} ProstT5 embeddings, colored by predicted substrate — "
            "explore the interactive version in the repo (5_umap/)",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=10.5, color="#555")
    order = [s for s in SUBSTRATE_COLORS if s in handles]
    ax.legend([handles[s] for s in order], [handles[s].get_label() for s in order],
              loc="upper right", fontsize=10, framealpha=0.9)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(OUT, dpi=150, bbox_inches="tight")
    plt.close()
    print("wrote", OUT)


if __name__ == "__main__":
    main()
