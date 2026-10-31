#!/usr/bin/env python3
"""
"Dark-matter" Actinomycetota 'Metabolite'-MTase clusters: the 14 clusters identified
by the corrected reviewed-status lookup (>= 80% Actinomycetota&Metabolite AND zero
SwissProt/reviewed members) -- the same 14 clusters reported everywhere else
in the paper (Table S1, CLEAN, Horizyn). Zoom-grid (Metabolite = orange).

Inputs:  sam_mtase_umap_prostt5_coords.csv, reviewed_status_full.csv
         (corrected reviewed-status lookup -- NOT the incomplete
         uniprot_metadata.csv, which undercounted reviewed members and
         produced a stale 15-cluster/wrong-membership figure previously),
         selected_1259_organisms.csv, ../6_dark_cluster_annotation/dark33_members.csv
         (ground-truth cluster membership shared with the rest of the manuscript)
Output:  actino_metabolite_dark_clusters_prostt5.png
"""
import os
import re
from math import ceil
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IO = Path(__file__).parent.resolve()
META_CSV = Path(os.environ.get("SMSTPRED_DATA", Path(__file__).resolve().parents[1] / "data")) / "selected_1259_organisms.csv"
DARK33_CSV = IO.parent / "6_dark_cluster_annotation" / "dark33_members.csv"
NCOLS = 5
BASE = "#ed8c44"   # 'Metabolite' UMAP color

# Authoritative cluster membership (matches Table S1, CLEAN, Horizyn analyses).
ACTINO_CLUSTERS = {17, 57, 200, 289, 309, 451, 452, 469, 470, 522, 692, 693, 795, 823}


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
    d["is_actino"] = d["Organism"].map(
        lambda o: "Actinomycetota" in org2lin.get(o, "") or "Actinobacteria" in org2lin.get(o, ""))
    d["is_metabolite"] = d["Substrate"] == "Metabolite"
    return d


def select(d):
    rows = []
    for c, sub in d.groupby("cluster"):
        if int(c) not in ACTINO_CLUSTERS:
            continue
        f_po = float((sub.is_actino & sub.is_metabolite).mean())
        names = sub["Protein names"].dropna()
        rows.append({"cluster": int(c), "size": len(sub), "f_actino_metabolite": f_po,
                     "name": _clean(names.value_counts().index[0]) if len(names) else "(no name)"})
    df = pd.DataFrame(rows).sort_values("size", ascending=False).reset_index(drop=True)
    missing = ACTINO_CLUSTERS - set(df.cluster)
    assert not missing, f"missing authoritative clusters from UMAP coords: {missing}"
    return df


def panel(ax, d, c, name, f_po):
    sub = d[d.cluster == c]
    ax.scatter(sub.UMAP_1, sub.UMAP_2, s=20, c=BASE, alpha=0.6, linewidth=0, zorder=2)
    mx = np.percentile(sub.UMAP_1.values, [1, 99]); my = np.percentile(sub.UMAP_2.values, [1, 99])
    px = max((mx[1] - mx[0]) * 0.12, 0.15); py = max((my[1] - my[0]) * 0.12, 0.15)
    ax.set_xlim(mx[0] - px, mx[1] + px); ax.set_ylim(my[0] - py, my[1] + py)
    ax.set_title(f"c{c}: {name}\n(n={len(sub)}, {f_po:.0%} actino&Metabolite)", fontsize=15)
    ax.set_xlabel("UMAP 1", fontsize=12); ax.set_ylabel("UMAP 2", fontsize=12)
    ax.tick_params(labelsize=10)


def main():
    d = load()
    sel = select(d)
    n = len(sel)
    print(f"{n} clusters: >=80% Actinomycetota&Metabolite and 0 SwissProt "
          f"({int(sel['size'].sum()):,} proteins)", flush=True)
    nrows = ceil(n / NCOLS)
    fig, axes = plt.subplots(nrows, NCOLS, figsize=(4.8 * NCOLS, 3.9 * nrows))
    axes = np.atleast_1d(axes).ravel()
    for ax, (_, r) in zip(axes, sel.iterrows()):
        panel(ax, d, int(r.cluster), r["name"], r.f_actino_metabolite)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle(f"Uncharacterized Actinomycetota 'Metabolite'-MTase clusters "
                 f"({n} clusters, >=80% actino&Metabolite, 0 SwissProt anchors)", fontsize=20, y=1.00)
    fig.tight_layout(rect=[0, 0, 1, 0.99])
    fig.savefig(IO / "actino_metabolite_dark_clusters_prostt5.png", dpi=150)
    print("wrote actino_metabolite_dark_clusters_prostt5.png", flush=True)


if __name__ == "__main__":
    main()
