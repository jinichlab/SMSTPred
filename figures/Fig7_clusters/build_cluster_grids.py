#!/usr/bin/env python3
"""Rebuild Figure 7 panels A-D (per-cluster UMAP scatter grids), stripped down:
no titles, no stats, no axis ticks/labels -- just the scatter and a big bold
cluster-ID number per subplot. No source script survived for the original
grids (only the final PNGs); descriptive per-cluster text (protein names,
SwissProt counts) turned out to drift against current UniProt data on
reconstruction, so it's dropped entirely rather than risk mismatched numbers.

Recipe:
  - Per-cluster membership: repo/5_umap/sam_mtase_umap_prostt5_coords.csv
    ("cluster" column), no extra lineage filter needed -- n matches the
    original figure exactly for every cluster checked.
  - Panel A/B (characterized, reviewed-anchored): SwissProt (reviewed) members
    highlighted gold via repo/5_umap/reviewed_status_full.csv.
  - Panel C/D (dark, uncharacterized) cluster lists/order:
    repo/6_dark_cluster_annotation/table2_dark_clusters.csv, sorted by n
    descending (matches the original reading order).
  - Badged clusters (Figure 7E numbered reactions): c551->1, c667->2, c110->3,
    c714->4 (plant dark set only).

Outputs (this dir): A_metazoa_protein_clusters.png, B_vertebrate_dna_clusters.png,
C_plant_other_dark_clusters.png, D_actino_other_dark_clusters.png
"""
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]      # repository root
COORDS = ROOT / "5_umap/sam_mtase_umap_prostt5_coords.csv"
REVIEWED = ROOT / "5_umap/reviewed_status_full.csv"
DARK_TABLE = ROOT / "6_dark_cluster_annotation/table2_dark_clusters.csv"
HERE = Path(__file__).resolve().parent

GOLD = "#FFD700"
BADGES = {551: "1", 667: "2", 110: "3", 714: "4"}


def load_data():
    coords = pd.read_csv(COORDS)
    rev = pd.read_csv(REVIEWED)
    return coords, rev


def bare_axes(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(1.2)


def draw_badge(ax, label):
    ax.add_patch(plt.Circle((0.08, 0.90), 0.1, transform=ax.transAxes,
                             facecolor="black", edgecolor="white", linewidth=1.4, zorder=10))
    ax.text(0.08, 0.90, label, transform=ax.transAxes, ha="center", va="center",
            color="white", fontsize=19, fontweight="bold", zorder=11)
    for s in ax.spines.values():
        s.set_linewidth(3.5)
        s.set_color("black")


def cluster_label(ax, cluster_id):
    ax.text(0.5, 1.02, str(cluster_id), transform=ax.transAxes, ha="center", va="bottom",
            fontsize=22, fontweight="bold")


def characterized_panel(ax, sub_members, rev, cluster_id, color):
    m = sub_members.merge(rev[["Entry", "Reviewed"]], on="Entry", how="left")
    is_rev = m["Reviewed"] == "reviewed"
    ax.scatter(m.loc[~is_rev, "UMAP_1"], m.loc[~is_rev, "UMAP_2"], s=26, c=color,
               alpha=0.55, linewidths=0, zorder=2)
    ax.scatter(m.loc[is_rev, "UMAP_1"], m.loc[is_rev, "UMAP_2"], s=110, c=GOLD,
               edgecolors="black", linewidths=0.9, zorder=4)
    cluster_label(ax, cluster_id)
    bare_axes(ax)


def dark_panel(ax, sub_members, cluster_id, color):
    ax.scatter(sub_members["UMAP_1"], sub_members["UMAP_2"], s=22, c=color,
               alpha=0.6, linewidths=0, zorder=2)
    cluster_label(ax, cluster_id)
    bare_axes(ax)
    if cluster_id in BADGES:
        draw_badge(ax, BADGES[cluster_id])


def main():
    coords, rev = load_data()

    def legend_handles(color):
        return [plt.Line2D([], [], marker="o", color="none", markerfacecolor=color,
                            alpha=0.55, markersize=13, label="cluster member (TrEMBL)"),
                plt.Line2D([], [], marker="o", color="none", markerfacecolor=GOLD,
                           markeredgecolor="black", markersize=15, label="SwissProt (reviewed)")]

    # ---- Panel A: Metazoa protein-methyltransferase clusters ----
    A_CLUSTERS = [399, 393, 441, 371, 114, 345]
    fig, axes = plt.subplots(2, 3, figsize=(15, 10.5), gridspec_kw={"wspace": 0.22, "hspace": 0.3})
    for cl, ax in zip(A_CLUSTERS, axes.flat):
        sub = coords[coords.cluster == cl]
        characterized_panel(ax, sub, rev, cl, "#238b45")
    fig.legend(handles=legend_handles("#238b45"), loc="lower center", ncol=2, fontsize=20,
               bbox_to_anchor=(0.5, -0.03), frameon=False)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(HERE / "A_metazoa_protein_clusters.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("[wrote] A_metazoa_protein_clusters.png")

    # ---- Panel B: Vertebrate DNA-methyltransferase clusters ----
    B_CLUSTERS = [158, 8, 64, 126]
    fig, axes = plt.subplots(2, 2, figsize=(11, 10.5), gridspec_kw={"wspace": 0.22, "hspace": 0.3})
    for cl, ax in zip(B_CLUSTERS, axes.flat):
        sub = coords[coords.cluster == cl]
        characterized_panel(ax, sub, rev, cl, "#6a51a3")
    fig.legend(handles=legend_handles("#6a51a3"), loc="lower center", ncol=2, fontsize=20,
               bbox_to_anchor=(0.5, -0.03), frameon=False)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(HERE / "B_vertebrate_dna_clusters.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("[wrote] B_vertebrate_dna_clusters.png")

    # ---- Panel C: plant dark clusters ----
    dark = pd.read_csv(DARK_TABLE)
    plant_ids = dark[dark.lineage == "Plant"].sort_values("n", ascending=False)["cluster"].tolist()
    fig, axes = plt.subplots(4, 5, figsize=(20, 15.5), gridspec_kw={"wspace": 0.3, "hspace": 0.35})
    for ax in axes.flat:
        ax.axis("off")
    for cl, ax in zip(plant_ids, axes.flat):
        ax.axis("on")
        sub = coords[coords.cluster == cl]
        dark_panel(ax, sub, cl, "#e6550d")
    fig.tight_layout()
    fig.savefig(HERE / "C_plant_other_dark_clusters.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("[wrote] C_plant_other_dark_clusters.png")

    # ---- Panel D: Actinomycetota dark clusters ----
    actino_ids = dark[dark.lineage == "Actinomycetota"].sort_values("n", ascending=False)["cluster"].tolist()
    fig, axes = plt.subplots(3, 5, figsize=(20, 11.5), gridspec_kw={"wspace": 0.3, "hspace": 0.35})
    for ax in axes.flat:
        ax.axis("off")
    for cl, ax in zip(actino_ids, axes.flat):
        ax.axis("on")
        sub = coords[coords.cluster == cl]
        dark_panel(ax, sub, cl, "#e6550d")
    fig.tight_layout()
    fig.savefig(HERE / "D_actino_other_dark_clusters.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("[wrote] D_actino_other_dark_clusters.png")


if __name__ == "__main__":
    main()
