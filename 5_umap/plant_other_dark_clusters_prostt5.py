#!/usr/bin/env python3
"""
"Dark-matter" plant 'Other'-MTase clusters: the 19 clusters identified by the
corrected reviewed-status lookup (>= 80% plant&Other AND zero SwissProt/
reviewed members) -- the same 19 clusters reported everywhere else in the
paper (Table 2/3, CLEAN, Horizyn). Zoom-grid (Other = orange), titled by the
dominant (TrEMBL) protein name.

Inputs:  sam_mtase_umap_prostt5_coords.csv, reviewed_status_full.csv
         (corrected reviewed-status lookup -- NOT the incomplete
         uniprot_metadata.csv, which undercounted reviewed members and
         produced a stale 20-cluster/wrong-membership figure previously),
         selected_1259_organisms.csv, ../6_dark_cluster_annotation/dark33_members.csv
         (ground-truth cluster membership shared with the rest of the manuscript)
Output:  plant_other_dark_clusters_prostt5.png
"""
import re
from math import ceil
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IO = Path(__file__).parent.resolve()
META_CSV = Path("/home/nroethler/Code/classifier/reference_files/selected_1259_organisms.csv")
DARK33_CSV = IO.parent / "6_dark_cluster_annotation" / "dark33_members.csv"
NCOLS = 5
BASE = "#e6550d"   # 'Other' UMAP color
BADGE_COLOR = "#1a1a1a"

# Authoritative cluster membership (matches Table 2/3, CLEAN, Horizyn analyses):
# the 19 dark33 clusters that are NOT Actinomycetota.
ACTINO_CLUSTERS = {17, 57, 200, 289, 309, 451, 452, 469, 470, 522, 692, 693, 795, 823}

# Clusters where CLEAN and Horizyn's dominant reaction family independently name
# the same (or closely related) specific chemistry -- cross-referenced to
# Figure 7E's numbered reaction list.
REACTION_BADGE = {551: "1", 667: "2", 110: "3", 714: "4"}


def _clean(name):
    name = re.sub(r"\s*\(EC [^)]*\)", "", str(name))
    return re.sub(r"\s*\((Fragment|isoform[^)]*)\)", "", name, flags=re.I)[:30]


def load():
    d = pd.read_csv(IO / "sam_mtase_umap_prostt5_coords.csv",
                    usecols=["Entry", "Organism", "Substrate", "UMAP_1", "UMAP_2", "cluster"])
    meta = pd.read_csv(IO / "reviewed_status_full.csv", dtype=str,
                       usecols=["Entry", "Reviewed", "Protein names"]).drop_duplicates("Entry")
    d = d.merge(meta, on="Entry", how="left")
    d["rev"] = d["Reviewed"] == "reviewed"
    org2lin = pd.read_csv(META_CSV).set_index("Organism")["Taxonomic lineage"].fillna("").to_dict()
    d["is_plant"] = d["Organism"].map(lambda o: "Viridiplantae" in org2lin.get(o, ""))
    d["is_other"] = d["Substrate"] == "Other"
    return d


def select(d):
    plant_clusters = sorted(set(pd.read_csv(DARK33_CSV).cluster.unique()) - ACTINO_CLUSTERS)
    rows = []
    for c, sub in d.groupby("cluster"):
        if int(c) not in plant_clusters:
            continue
        names = sub["Protein names"].dropna()
        rows.append({"cluster": int(c), "size": len(sub),
                     "f_plant_other": float((sub.is_plant & sub.is_other).mean()),
                     "name": _clean(names.value_counts().index[0]) if len(names) else "(no name)"})
    df = pd.DataFrame(rows).sort_values("size", ascending=False).reset_index(drop=True)
    missing = set(plant_clusters) - set(df.cluster)
    assert not missing, f"missing authoritative clusters from UMAP coords: {missing}"
    return df


def panel(ax, d, c, name, f_po):
    sub = d[d.cluster == c]
    ax.scatter(sub.UMAP_1, sub.UMAP_2, s=20, c=BASE, alpha=0.6, linewidth=0, zorder=2)
    mx = np.percentile(sub.UMAP_1.values, [1, 99]); my = np.percentile(sub.UMAP_2.values, [1, 99])
    px = max((mx[1] - mx[0]) * 0.12, 0.15); py = max((my[1] - my[0]) * 0.12, 0.15)
    ax.set_xlim(mx[0] - px, mx[1] + px); ax.set_ylim(my[0] - py, my[1] + py)
    ax.set_title(f"c{c}: {name}\n(n={len(sub)}, {f_po:.0%} plant&Other)", fontsize=15)
    ax.set_xlabel("UMAP 1", fontsize=12); ax.set_ylabel("UMAP 2", fontsize=12)
    ax.tick_params(labelsize=10)

    badge = REACTION_BADGE.get(c)
    if badge:
        for spine in ax.spines.values():
            spine.set_edgecolor(BADGE_COLOR); spine.set_linewidth(3)
        ax.text(0.14, 0.86, badge, transform=ax.transAxes, ha="center", va="center",
                fontsize=22, fontweight="bold", color="white", zorder=6,
                bbox=dict(boxstyle="circle,pad=0.35", fc=BADGE_COLOR, ec="white", lw=1.2))


def main():
    d = load()
    sel = select(d)
    n = len(sel)
    print(f"{n} clusters: 0 SwissProt, plant&Other "
          f"{sel.f_plant_other.min():.0%}-{sel.f_plant_other.max():.0%}", flush=True)
    nrows = ceil(n / NCOLS)
    fig, axes = plt.subplots(nrows, NCOLS, figsize=(4.8 * NCOLS, 3.9 * nrows))
    axes = np.atleast_1d(axes).ravel()
    for ax, (_, r) in zip(axes, sel.iterrows()):
        panel(ax, d, int(r.cluster), r["name"], r.f_plant_other)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle(f"Uncharacterized plant 'Other'-MTase clusters "
                 f"({n} clusters, {sel.f_plant_other.min():.0%}+ plant&Other, 0 SwissProt anchors)",
                 fontsize=20, y=1.00)
    fig.tight_layout(rect=[0, 0, 1, 0.99])
    fig.savefig(IO / "plant_other_dark_clusters_prostt5.png", dpi=150)
    print("wrote plant_other_dark_clusters_prostt5.png", flush=True)


if __name__ == "__main__":
    main()
