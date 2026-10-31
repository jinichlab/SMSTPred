#!/usr/bin/env python3
"""Reproduce the dark-cluster definition used in the manuscript.

A dark cluster is an HDBSCAN cluster of the ProstT5 UMAP that contains no reviewed (Swiss-Prot) member and in which at least 80 percent
of the proteins are BOTH from one lineage (Actinomycetota or plants = Viridiplantae) AND predicted 'Metabolite' (single-class call). The
80 percent applies to the joint condition within each cluster, not to the lineage and the class separately.

Inputs: sam_mtase_umap_prostt5_coords.csv (cluster, predicted substrate), reviewed_status_full.csv, data/selected_1259_organisms.csv (lineage).
Outputs (next to this script): dark_cluster_definition.csv (every cluster with no reviewed member and its lineage / 'Metabolite' shares) and
dark_cluster_definition.txt (the counts quoted in the Methods and a check against 6_dark_cluster_annotation/dark33_members.csv).
"""
import os
import re
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
DATA = Path(os.environ.get("SMSTPRED_DATA", REPO / "data"))
THRESHOLD = 0.80
LINEAGES = {"Actinomycetota": "Actinomycetota", "Plant": "Viridiplantae"}


def canon(name):
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9_\-]+", "", re.sub(r"\s+", "_", str(name).strip()))).strip("_")


coords = pd.read_csv(HERE / "sam_mtase_umap_prostt5_coords.csv", low_memory=False)
rev = pd.read_csv(HERE / "reviewed_status_full.csv", usecols=["Entry", "Reviewed"]).drop_duplicates("Entry")
meta = pd.read_csv(DATA / "selected_1259_organisms.csv")
meta["key"] = meta.Organism.map(canon)
coords["key"] = coords.Organism.map(canon)
d = coords.merge(rev, on="Entry", how="left").merge(meta[["key", "Taxonomic lineage"]].drop_duplicates("key"), on="key", how="left")
d["reviewed"] = d.Reviewed.fillna("").astype(str).str.lower().eq("reviewed")
d["is_metabolite"] = d.Substrate.eq("Metabolite")
for name, tax in LINEAGES.items():
    d[f"in_{name}"] = d["Taxonomic lineage"].fillna("").str.contains(tax)
d = d[d.cluster >= 0]                                    # -1 is HDBSCAN noise

rows = []
for c, g in d.groupby("cluster"):
    r = dict(cluster=int(c), n=len(g), n_reviewed=int(g.reviewed.sum()), share_metabolite=g.is_metabolite.mean())
    for name in LINEAGES:
        r[f"share_{name}"] = g[f"in_{name}"].mean()
        r[f"share_{name}_and_metabolite"] = (g[f"in_{name}"] & g.is_metabolite).mean()
    rows.append(r)
cl = pd.DataFrame(rows)
no_rev = cl[cl.n_reviewed == 0].copy()
no_rev["lineage"] = ""
for name in LINEAGES:
    no_rev.loc[no_rev[f"share_{name}_and_metabolite"] >= THRESHOLD, "lineage"] = name
dark = no_rev[no_rev.lineage != ""]
no_rev.sort_values(["lineage", "n"], ascending=[False, False]).to_csv(HERE / "dark_cluster_definition.csv", index=False, float_format="%.4f")

members = pd.read_csv(REPO / "6_dark_cluster_annotation" / "dark33_members.csv")
same = set(dark.cluster) == set(members.cluster)
n_dark = int(dark.n.sum())
lines = [f"HDBSCAN clusters (noise excluded): {len(cl)}",
         f"clusters with no reviewed member: {len(no_rev)}",
         f"dark clusters (>= {THRESHOLD:.0%} of members in one lineage AND predicted Metabolite): {len(dark)} "
         f"({(dark.lineage == 'Actinomycetota').sum()} Actinomycetota, {(dark.lineage == 'Plant').sum()} plant)",
         f"proteins in the dark clusters: {n_dark:,} of {int(cl.n.sum()):,} clustered proteins ({n_dark / cl.n.sum():.1%})",
         f"same 33 clusters as 6_dark_cluster_annotation/dark33_members.csv: {same}",
         f"(for comparison, requiring the lineage share and the 'Metabolite' share to be >= {THRESHOLD:.0%} separately instead of jointly "
         f"gives {int(((no_rev.share_metabolite >= THRESHOLD) & ((no_rev.share_Actinomycetota >= THRESHOLD) | (no_rev.share_Plant >= THRESHOLD))).sum())} clusters)"]
(HERE / "dark_cluster_definition.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
