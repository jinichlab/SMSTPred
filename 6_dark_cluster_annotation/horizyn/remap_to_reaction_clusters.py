#!/usr/bin/env python3
"""Remap each dark protein-cluster's per-member top-1 vote from a raw Rhea
reaction ID onto its REACTION-CLUSTER label (from cluster_rhea_reactions.py),
then recompute the dominant-reaction / diversity summary at the family level
instead of the raw-ID level. Also computes each reaction cluster's own
orphan-ness (fraction of its member reactions with no annotated enzyme).
"""
import os
import numpy as np
import pandas as pd

HERE = os.environ.get("HORIZYN_DIR", os.path.dirname(os.path.abspath(__file__)))
TOPK_CSV = f"{HERE}/horizyn_full_topk.csv"
RXN_CLUSTERS_CSV = f"{HERE}/rhea_reaction_clusters.csv"
OUT_SUMMARY = f"{HERE}/protein_clusters_by_reaction_family.csv"
OUT_RXN_FAMILIES = f"{HERE}/reaction_family_summary.csv"


def tier(h):
    if h < 0.15:
        return "Monofunctional"
    if h < 0.30:
        return "Dominant chemistry"
    if h < 0.45:
        return "Mixed"
    return "Heterogeneous"


def main():
    top = pd.read_csv(TOPK_CSV)
    r1 = top[top["rank"] == 1].copy()

    rxn = pd.read_csv(RXN_CLUSTERS_CSV)  # rhea_lr_id, reaction_cluster, is_sam_methylation, has_enzyme, reaction_smiles
    rxn_lookup = rxn.set_index("rhea_lr_id")

    r1 = r1.merge(rxn[["rhea_lr_id", "reaction_cluster"]], on="rhea_lr_id", how="left")
    n_unmapped = r1.reaction_cluster.isna().sum()
    print(f"[remap] {n_unmapped}/{len(r1)} rank-1 votes have no reaction-cluster label "
          f"(HDBSCAN noise, label -1 counts as mapped; NaN = missing from rxn table)")

    # --- per reaction-family (cluster) summary: size, %orphan, sample equations ---
    fam_rows = []
    for fc, g in rxn.groupby("reaction_cluster"):
        n = len(g)
        n_orphan = int((~g.has_enzyme).sum())
        n_sam = int(g.is_sam_methylation.sum())
        fam_rows.append(dict(reaction_cluster=fc, n_reactions=n, n_orphan=n_orphan,
                              pct_orphan=round(n_orphan/n*100, 1), n_sam_methylation=n_sam,
                              sample_reactions="; ".join(g.reaction_smiles.head(3))))
    fam_df = pd.DataFrame(fam_rows).sort_values("n_reactions", ascending=False)
    fam_df.to_csv(OUT_RXN_FAMILIES, index=False)
    print(f"[remap] {len(fam_df)} reaction families (incl. -1 noise) -> {OUT_RXN_FAMILIES}")

    # --- per protein-cluster: dominant reaction FAMILY, entropy over families ---
    rows = []
    for pcluster, g in r1.groupby("cluster"):
        n = len(g)
        g_valid = g.dropna(subset=["reaction_cluster"])
        counts = g_valid.reaction_cluster.value_counts()
        if len(counts) == 0:
            continue
        p = counts / len(g_valid)
        H = max(-(p * np.log2(p)).sum() / np.log2(len(g_valid)), 0.0) if len(g_valid) > 1 else 0.0
        top_fam = counts.index[0]
        mode_share = counts.iloc[0] / n
        fam_info = fam_df[fam_df.reaction_cluster == top_fam].iloc[0]
        rows.append(dict(
            cluster=pcluster, n=n, tier=tier(H), entropy=round(H, 3),
            n_distinct_families=g_valid.reaction_cluster.nunique(),
            dominant_family=int(top_fam), dominant_family_share=round(mode_share, 3),
            dominant_family_size=fam_info.n_reactions, dominant_family_pct_orphan=fam_info.pct_orphan,
            dominant_family_sample=fam_info.sample_reactions,
        ))
    out = pd.DataFrame(rows).sort_values("entropy")
    out.to_csv(OUT_SUMMARY, index=False)
    print(f"[remap] wrote {OUT_SUMMARY}")
    print(out[["cluster","n","tier","entropy","n_distinct_families",
               "dominant_family","dominant_family_share","dominant_family_pct_orphan"]].to_string(index=False))


if __name__ == "__main__":
    main()
