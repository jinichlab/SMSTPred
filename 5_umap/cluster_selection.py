#!/usr/bin/env python3
"""Which clusters of a UMAP are the dark clusters, the metazoan protein-substrate clusters and the vertebrate DNA-substrate clusters.

  python cluster_selection.py [--release <folder>] [--embedding prott5]     (default release: $RELEASE_DIR, else the newest release_* of $UNIPROT_PULL_DIR)

Reads <release>/umap/sam_mtase_umap_<embedding>_coords.csv (global_umap.py), the reviewed (Swiss-Prot) status and protein names of <release>/reference_proteomes/*.csv and
the lineage flags of data/selected_1259_organisms.csv. The rules (one HDBSCAN cluster each; "call" = the single class a protein is called):
  dark cluster                 no reviewed member AND >= 80 % of the members are BOTH in one lineage (Actinomycetota, or plants = Viridiplantae) AND called Metabolite
                               (the definition of 5_umap/dark_cluster_definition.py: the 80 % applies to the joint condition)
  metazoan protein cluster     >= 90 % of the members metazoan AND >= 90 % called Protein, and >= 1 reviewed member (an anchor)
  vertebrate DNA cluster       >= 50 % of the members called DNA AND >= 30 vertebrate members; OR (the anchor clause) >= 30 % called DNA AND >= 20 vertebrates AND a reviewed
                               vertebrate member whose name has 'DNA' and 'methyltransferase' (catches the Dnmt3a cluster of the ProtT5 map, which is only 43 % DNA-called)
The last two rules were written to reproduce the clusters of the archived Fig 7A/B (6 and 4 clusters of the ProstT5 map; see test_cluster_selection.py), no rule for them
survived; on that map the anchor clause adds nothing. The 30 vertebrates sit between the archived clusters (31 and 34) and the next candidates (29 and 24), so that cut is a fit:
the near misses are listed in the report.
Writes to <release>/umap/: cluster_table.csv (every cluster), dark_clusters.csv, metazoa_protein_clusters.csv, vertebrate_dna_clusters.csv,
cluster_crosswalk_vs_prostt5.csv (the archived clusters of Fig 7/8 -> the new cluster that shares the most proteins) and cluster_selection_report.txt.
"""
import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
META_CSV = Path(os.environ.get("SMSTPRED_DATA", REPO / "data")) / "selected_1259_organisms.csv"
DARK_THRESHOLD = 0.80
LINEAGES = {"Actinomycetota": "Actinomycetota", "Plant": "Viridiplantae"}
METAZOA_PROTEIN = dict(metazoa=0.90, protein=0.90, min_reviewed=1)
VERTEBRATE_DNA = dict(dna=0.50, min_vertebrates=30, anchor_dna=0.30, anchor_vertebrates=20)
# the archived (ProstT5-map) clusters of Fig 7A/B/C/D and the four of Fig 8, for the crosswalk
ARCHIVED_A = [399, 393, 441, 371, 114, 345]
ARCHIVED_B = [158, 8, 64, 126]
ARCHIVED_FIG8 = [551, 667, 110, 714]


def clean_name(n):
    return re.sub(r"\s*\(EC [^)]*\)", "", str(n)).strip()[:60]


def annotate(coords, reviewed, meta):
    """coords: Entry, Organism, Substrate, cluster; reviewed: Entry, Reviewed, Protein names; meta: Organism, Taxonomic lineage, is_metazoa, is_vertebrate."""
    d = coords.merge(reviewed.drop_duplicates("Entry"), on="Entry", how="left").merge(
        meta.drop_duplicates("Organism")[["Organism", "Taxonomic lineage", "is_metazoa", "is_vertebrate"]], on="Organism", how="left")
    lin = d["Taxonomic lineage"].fillna("")
    d["reviewed"] = d["Reviewed"].fillna("").astype(str).str.lower().eq("reviewed")
    d["metazoa"] = d["is_metazoa"].fillna(False).astype(bool)
    d["vertebrate"] = d["is_vertebrate"].fillna(False).astype(bool)
    nm = d["Protein names"].fillna("").astype(str).str.lower()
    d["rev_vert_dnmt"] = d.reviewed & d.vertebrate & nm.str.contains("dna") & nm.str.contains("methyltransferase")
    for name, tax in LINEAGES.items():
        d[f"in_{name}"] = lin.str.contains(tax)
    for lab in ["Metabolite", "Protein", "RNA", "DNA"]:
        d[f"is_{lab}"] = d["Substrate"].eq(lab)
    return d[d.cluster >= 0]


def cluster_table(d):
    g = d.groupby("cluster")
    t = pd.DataFrame({"n": g.size(), "n_reviewed": g.reviewed.sum().astype(int), "n_vertebrate": g.vertebrate.sum().astype(int), "n_rev_vert_dnmt": g.rev_vert_dnmt.sum().astype(int),
                      "share_metazoa": g.metazoa.mean(), "share_vertebrate": g.vertebrate.mean()})
    for lab in ["Metabolite", "Protein", "RNA", "DNA"]:
        t[f"share_{lab}"] = g[f"is_{lab}"].mean()
    for name in LINEAGES:
        t[f"share_{name}"] = g[f"in_{name}"].mean()
        t[f"share_{name}_and_Metabolite"] = d.assign(x=d[f"in_{name}"] & d.is_Metabolite).groupby("cluster").x.mean()
    names = d[d.reviewed].dropna(subset=["Protein names"]).assign(nm=lambda x: x["Protein names"].map(clean_name)).groupby("cluster").nm.agg(
        lambda s: "; ".join(f"{k} (x{v})" if v > 1 else k for k, v in s.value_counts().head(3).items()))
    t["reviewed_names"] = names
    t["reviewed_names"] = t["reviewed_names"].fillna("")
    return t.reset_index()


def select_dark(t):
    nr = t[t.n_reviewed == 0].copy()
    nr["lineage"] = ""
    for name in LINEAGES:
        nr.loc[nr[f"share_{name}_and_Metabolite"] >= DARK_THRESHOLD, "lineage"] = name
    return nr[nr.lineage != ""].sort_values(["lineage", "n"], ascending=[True, False])


def select_metazoa_protein(t):
    p = METAZOA_PROTEIN
    return t[(t.share_metazoa >= p["metazoa"]) & (t.share_Protein >= p["protein"]) & (t.n_reviewed >= p["min_reviewed"])].sort_values("n", ascending=False)


def select_vertebrate_dna(t):
    p = VERTEBRATE_DNA
    base = (t.share_DNA >= p["dna"]) & (t.n_vertebrate >= p["min_vertebrates"])
    anchor = (t.share_DNA >= p["anchor_dna"]) & (t.n_vertebrate >= p["anchor_vertebrates"]) & (t.n_rev_vert_dnmt >= 1)
    return t[base | anchor].sort_values("n_vertebrate", ascending=False)


def crosswalk(new_coords, old_coords, sets):
    """For each archived cluster in sets {name: [ids]}: the new cluster sharing the most proteins. jaccard = shared / union; old_in_new = shared / archived size."""
    new = new_coords[new_coords.cluster >= 0].groupby("cluster").Entry.agg(set)
    old = old_coords[old_coords.cluster >= 0].groupby("cluster").Entry.agg(set)
    inv = new_coords[new_coords.cluster >= 0].set_index("Entry").cluster
    rows = []
    for setname, ids in sets.items():
        for c in ids:
            if c not in old.index:
                continue
            members = old[c]
            hits = inv.reindex(list(members)).dropna().astype(int).value_counts()
            if hits.empty:
                rows.append({"set": setname, "old_cluster": c, "old_n": len(members)})
                continue
            nc, shared = int(hits.index[0]), int(hits.iloc[0])
            rows.append({"set": setname, "old_cluster": c, "old_n": len(members), "new_cluster": nc, "new_n": len(new[nc]), "shared": shared,
                         "jaccard": shared / len(members | new[nc]), "old_in_new": shared / len(members)})
    return pd.DataFrame(rows)


def load_release(rel, emb):
    coords = pd.read_csv(rel / "umap" / f"sam_mtase_umap_{emb}_coords.csv", low_memory=False)
    keep = set(coords.Entry)
    parts = []
    slim = rel / "embedding_inputs" / "protein_annotations.csv.gz"      # the archived record ships this slim table instead of the 3 GB of proteome tables
    if (rel / "reference_proteomes").is_dir() and any((rel / "reference_proteomes").glob("*.csv")):
        for f in sorted((rel / "reference_proteomes").glob("*.csv")):
            x = pd.read_csv(f, usecols=["Entry", "Reviewed", "Protein names"], low_memory=False)
            parts.append(x[x.Entry.isin(keep)])
    else:
        x = pd.read_csv(slim, usecols=["Entry", "Reviewed", "Protein names"], low_memory=False)
        parts.append(x[x.Entry.isin(keep)])
    return coords, pd.concat(parts, ignore_index=True), pd.read_csv(META_CSV)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release")
    ap.add_argument("--embedding", default="prott5")
    a = ap.parse_args()
    sys.path.insert(0, str(REPO / "1_clustering"))
    from mmseqs_clusters import newest_release
    rel = Path(a.release or os.environ.get("RELEASE_DIR") or newest_release(Path(os.environ.get("UNIPROT_PULL_DIR", REPO / "data" / "uniprot_pull")), os.environ.get("RELEASE")))
    out = rel / "umap"
    coords, reviewed, meta = load_release(rel, a.embedding)
    d = annotate(coords, reviewed, meta)
    t = cluster_table(d)
    dark, mp, vd = select_dark(t), select_metazoa_protein(t), select_vertebrate_dna(t)
    t.to_csv(out / "cluster_table.csv", index=False, float_format="%.4f")
    dark.to_csv(out / "dark_clusters.csv", index=False, float_format="%.4f")
    mp.to_csv(out / "metazoa_protein_clusters.csv", index=False, float_format="%.4f")
    vd.to_csv(out / "vertebrate_dna_clusters.csv", index=False, float_format="%.4f")

    L = [f"{len(d):,} proteins in {len(t)} clusters ({coords.shape[0] - len(d):,} unclustered), {int(t.n_reviewed.sum()):,} reviewed",
         f"clusters with no reviewed member: {int((t.n_reviewed == 0).sum())}",
         f"DARK clusters (>= {DARK_THRESHOLD:.0%} of members one lineage AND called Metabolite, 0 reviewed): {len(dark)} "
         f"({int((dark.lineage == 'Actinomycetota').sum())} Actinomycetota, {int((dark.lineage == 'Plant').sum())} plant); {int(dark.n.sum()):,} proteins = {dark.n.sum() / t.n.sum():.1%} of the clustered ones",
         f"METAZOAN PROTEIN clusters ({METAZOA_PROTEIN}): {len(mp)}",
         f"VERTEBRATE DNA clusters ({VERTEBRATE_DNA}): {len(vd)}", ""]
    L.append("metazoan protein clusters (cluster, n, reviewed, metazoan share, Protein share, reviewed names):")
    for _, r in mp.iterrows():
        L.append(f"  c{int(r.cluster)}: n={int(r.n)}, reviewed {int(r.n_reviewed)}, metazoa {r.share_metazoa:.0%}, Protein {r.share_Protein:.0%}  {r.reviewed_names}")
    near = t[(t.share_metazoa >= .8) & (t.share_Protein >= .8) & (t.n_reviewed >= 1) & ~t.cluster.isin(mp.cluster)]
    L.append(f"  near misses (>= 80 % on both): {[(int(r.cluster), f'{r.share_metazoa:.0%}/{r.share_Protein:.0%}') for _, r in near.iterrows()]}")
    L.append("vertebrate DNA clusters (cluster, n, vertebrates, DNA share, reviewed names):")
    for _, r in vd.iterrows():
        L.append(f"  c{int(r.cluster)}: n={int(r.n)}, vertebrates {int(r.n_vertebrate)}, DNA {r.share_DNA:.0%}, reviewed {int(r.n_reviewed)}  {r.reviewed_names}")
    near = t[(t.share_DNA >= .3) & (t.n_vertebrate >= 15) & ~t.cluster.isin(vd.cluster)]
    L.append(f"  near misses (DNA >= 30 %, >= 15 vertebrates, not selected): {[(int(r.cluster), int(r.n_vertebrate)) for _, r in near.iterrows()]}")

    old_path, dark33 = HERE / "sam_mtase_umap_prostt5_coords.csv", REPO / "6_dark_cluster_annotation" / "table2_dark_clusters.csv"
    if old_path.exists() and dark33.exists():
        old = pd.read_csv(old_path, usecols=["Entry", "cluster"])
        t2 = pd.read_csv(dark33)
        sets = {"A_metazoa_protein": ARCHIVED_A, "B_vertebrate_dna": ARCHIVED_B, "dark_Actinomycetota": t2[t2.lineage == "Actinomycetota"].cluster.tolist(),
                "dark_Plant": t2[t2.lineage == "Plant"].cluster.tolist(), "Fig8": ARCHIVED_FIG8}
        cw = crosswalk(coords, old, sets)
        selected = {"A_metazoa_protein": set(mp.cluster), "B_vertebrate_dna": set(vd.cluster), "dark_Actinomycetota": set(dark[dark.lineage == "Actinomycetota"].cluster),
                    "dark_Plant": set(dark[dark.lineage == "Plant"].cluster), "Fig8": set(dark.cluster)}
        cw["new_cluster_selected"] = [int(r.new_cluster) in selected[r.set] if pd.notna(r.get("new_cluster")) else False for _, r in cw.iterrows()]
        cw.to_csv(out / "cluster_crosswalk_vs_prostt5.csv", index=False, float_format="%.3f")
        L += ["", "crosswalk to the archived ProstT5 map (an archived cluster counts as found when one new cluster has Jaccard >= 0.5 with it AND passes the same rule):"]
        for name, g in cw.groupby("set", sort=False):
            found = g[(g.jaccard >= 0.5) & g.new_cluster_selected]
            L.append(f"  {name}: {len(found)} of {len(g)} found; median Jaccard {g.jaccard.median():.2f}, median share of the old members in the best new cluster {g.old_in_new.median():.0%}")
        old_dark_found = cw[cw.set.str.startswith("dark") & (cw.jaccard >= 0.5) & cw.new_cluster_selected].new_cluster.nunique()
        L.append(f"  new dark clusters matched by an archived dark cluster: {old_dark_found} of {len(dark)}")
    (out / "cluster_selection_report.txt").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
