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
  - Outlined clusters (the four examined in Figure 8): c551, c667, c110, c714
    (plant dark set only). Thick outline, no number.

Outputs (this dir): A_metazoa_protein_clusters.png, B_vertebrate_dna_clusters.png,
C_plant_metabolite_dark_clusters.png, D_actino_metabolite_dark_clusters.png

  python build_cluster_grids.py                       the archived ProstT5 map and its hard-coded cluster lists; writes here
  python build_cluster_grids.py --release [<folder>]  a refresh release's ProtT5 map: the clusters are the ones selected by 5_umap/cluster_selection.py
                                                      (<release>/umap/{metazoa_protein,vertebrate_dna,dark}_clusters.csv, reviewed status from
                                                      <release>/reference_proteomes); writes <release>/umap/fig7/. The panels of the four clusters of Fig 8 are outlined
                                                      where the archived cluster has a counterpart (Jaccard >= 0.5, cluster_crosswalk_vs_prostt5.csv).
"""
import argparse
import os
import sys
from math import ceil
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

GOLD = "#f4d23b"
OUTLINED = {551, 667, 110, 714}      # the four dark clusters examined in Figure 8


def load_data():
    coords = pd.read_csv(COORDS)
    rev = pd.read_csv(REVIEWED)
    return coords, rev


def bare_axes(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(1.2)


def outline_panel(ax):
    for s_ in ax.spines.values():
        s_.set_linewidth(3.5)
        s_.set_color("black")


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
    if cluster_id in OUTLINED:
        outline_panel(ax)


def grid_axes(n, ncols, cell=(4.0, 3.5)):
    nrows = ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(cell[0] * ncols, cell[1] * nrows), gridspec_kw={"wspace": 0.3, "hspace": 0.35})
    axes = axes.ravel() if hasattr(axes, "ravel") else [axes]
    for ax in axes:
        ax.axis("off")
    return fig, axes


def main_release(release_arg):
    global OUTLINED
    sys.path.insert(0, str(ROOT / "1_clustering"))
    from mmseqs_clusters import newest_release
    rel = Path(release_arg or os.environ.get("RELEASE_DIR") or newest_release(Path(os.environ.get("UNIPROT_PULL_DIR", ROOT / "data" / "uniprot_pull")), os.environ.get("RELEASE")))
    u = rel / "umap"
    out = u / "fig7"
    out.mkdir(exist_ok=True)
    coords = pd.read_csv(u / "sam_mtase_umap_prott5_coords.csv", low_memory=False)
    keep = set(coords.Entry)
    files = sorted((rel / "reference_proteomes").glob("*.csv")) if (rel / "reference_proteomes").is_dir() else []
    if files:
        rev = pd.concat([pd.read_csv(f, usecols=["Entry", "Reviewed"], low_memory=False) for f in files]).drop_duplicates("Entry")
    else:      # the archived record ships a slim table instead of the per-proteome tables
        rev = pd.read_csv(rel / "embedding_inputs" / "protein_annotations.csv.gz", usecols=["Entry", "Reviewed"], low_memory=False).drop_duplicates("Entry")
    rev = rev[rev.Entry.isin(keep)]
    mp = pd.read_csv(u / "metazoa_protein_clusters.csv").sort_values("n", ascending=False).cluster.tolist()
    vd = pd.read_csv(u / "vertebrate_dna_clusters.csv").sort_values("n", ascending=False).cluster.tolist()
    dark = pd.read_csv(u / "dark_clusters.csv")
    plant = dark[dark.lineage == "Plant"].sort_values("n", ascending=False).cluster.tolist()
    actino = dark[dark.lineage == "Actinomycetota"].sort_values("n", ascending=False).cluster.tolist()
    OUTLINED = set()
    cwp = u / "cluster_crosswalk_vs_prostt5.csv"
    f8p = rel / "fig8" / "fig8_clusters.csv"       # the clusters Fig 8 shows (written by make_fig8_reaction_profiles.py --release): they are outlined here
    if f8p.exists():
        OUTLINED = set(pd.read_csv(f8p).cluster.astype(int))
        print(f"Fig 8 shows {sorted(OUTLINED)} -> outlined where they are in the dark Plant set")
    elif cwp.exists():
        cw = pd.read_csv(cwp)
        f8 = cw[(cw.set == "Fig8") & (cw.jaccard >= 0.5) & cw.new_cluster_selected]
        OUTLINED = set(f8.new_cluster.astype(int))
        print(f"Fig 8 clusters with a counterpart on the new map: {len(f8)} of 4 -> outlined {sorted(OUTLINED)}")

    def legend_handles(color):
        return [plt.Line2D([], [], marker="o", color="none", markerfacecolor=color, alpha=0.55, markersize=13, label="cluster member (TrEMBL)"),
                plt.Line2D([], [], marker="o", color="none", markerfacecolor=GOLD, markeredgecolor="black", markersize=15, label="SwissProt (reviewed)")]

    for ids, color, ncols, name in ((mp, "#618d43", 4, "A_metazoa_protein_clusters"), (vd, "#6a51a3", 3, "B_vertebrate_dna_clusters")):
        fig, axes = grid_axes(len(ids), ncols, (4.6, 4.2))
        for cl, ax in zip(ids, axes):
            ax.axis("on")
            characterized_panel(ax, coords[coords.cluster == cl], rev, cl, color)
        fig.legend(handles=legend_handles(color), loc="lower center", ncol=2, fontsize=20, bbox_to_anchor=(0.5, -0.03), frameon=False)
        fig.tight_layout(rect=[0, 0.04, 1, 1])
        fig.savefig(out / f"{name}.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"[wrote] {name}.png ({len(ids)} clusters)")
    for ids, ncols, name in ((plant, 6, "C_plant_metabolite_dark_clusters"), (actino, 6, "D_actino_metabolite_dark_clusters")):
        fig, axes = grid_axes(len(ids), ncols, (3.9, 3.5))
        for cl, ax in zip(ids, axes):
            ax.axis("on")
            dark_panel(ax, coords[coords.cluster == cl], cl, "#ed8c44")
        fig.tight_layout()
        fig.savefig(out / f"{name}.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"[wrote] {name}.png ({len(ids)} clusters)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", nargs="?", const="", default=None)
    args = ap.parse_args()
    if args.release is not None:
        return main_release(args.release)
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
        characterized_panel(ax, sub, rev, cl, "#618d43")
    fig.legend(handles=legend_handles("#618d43"), loc="lower center", ncol=2, fontsize=20,
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
        dark_panel(ax, sub, cl, "#ed8c44")
    fig.tight_layout()
    fig.savefig(HERE / "C_plant_metabolite_dark_clusters.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("[wrote] C_plant_metabolite_dark_clusters.png")

    # ---- Panel D: Actinomycetota dark clusters ----
    actino_ids = dark[dark.lineage == "Actinomycetota"].sort_values("n", ascending=False)["cluster"].tolist()
    fig, axes = plt.subplots(3, 5, figsize=(20, 11.5), gridspec_kw={"wspace": 0.3, "hspace": 0.35})
    for ax in axes.flat:
        ax.axis("off")
    for cl, ax in zip(actino_ids, axes.flat):
        ax.axis("on")
        sub = coords[coords.cluster == cl]
        dark_panel(ax, sub, cl, "#ed8c44")
    fig.tight_layout()
    fig.savefig(HERE / "D_actino_metabolite_dark_clusters.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("[wrote] D_actino_metabolite_dark_clusters.png")


if __name__ == "__main__":
    main()
