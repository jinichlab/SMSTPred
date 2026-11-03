#!/usr/bin/env python3
"""Per dark protein-cluster Horizyn summary: the reaction family its members vote for, and the mean score of their best match.

  python final_family_table.py --topk <horizyn_full_topk.csv> --out <folder> [--rxn-clusters rhea_reaction_clusters.csv]

Each protein votes with its rank-1 reaction (horizyn_full_topk.csv, run_horizyn_predict_full.py), which is mapped to its reaction family (rhea_reaction_clusters.csv,
cluster_rhea_reactions.py; -1 = HDBSCAN noise, which is NOT a family). Per protein cluster:
  n                          proteins in the cluster
  pct_noise_votes            percent of the members whose rank-1 reaction is noise (-1)
  coherent_dominant_family   the family with the most votes among the non-noise votes
  coherent_family_share      its votes as a share of ALL members (noise votes included, so a cluster with 76 % noise votes cannot have a share above 24 %)
  coherent_family_size / _pct_orphan / _sample   the family's reactions, the percent of them without an annotated enzyme, the first three reaction SMILES (cut to 200 characters, as in the archived table)
  candidate_score            share x orphan fraction (a high value = most members vote for a family of reactions that no enzyme is known for)
and avg_top1_score = the mean cosine score of the members' own rank-1 match. Writes protein_clusters_family_final.csv and avg_top1_score.csv to --out.
Reproduces the archived tables of the 33 dark clusters from the archived top-50 file (test_final_family_table.py).
"""
import argparse
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent


def family_table(rxn):
    rows = []
    for fc, g in rxn.groupby("reaction_cluster"):
        n = len(g)
        rows.append(dict(reaction_cluster=fc, n_reactions=n, pct_orphan=round(int((~g.has_enzyme).sum()) / n * 100, 1), sample_reactions="; ".join(g.reaction_smiles.head(3))))
    return pd.DataFrame(rows).set_index("reaction_cluster")


def summarize(topk, rxn):
    fam = family_table(rxn)
    r1 = topk[topk["rank"] == 1].merge(rxn[["rhea_lr_id", "reaction_cluster"]], on="rhea_lr_id", how="left")
    rows, avg = [], []
    for c, g in r1.groupby("cluster"):
        n = len(g)
        votes = g.reaction_cluster.dropna().astype(int).value_counts()
        coherent = votes[votes.index != -1]
        noise = int(votes.get(-1, 0))
        avg.append({"cluster": c, "avg_top1_score": round(g.score.mean(), 3)})
        if coherent.empty:
            continue
        f = int(coherent.index[0])
        share = round(coherent.iloc[0] / n, 3)
        rows.append(dict(cluster=c, n=n, pct_noise_votes=round(noise / n * 100, 1), coherent_dominant_family=f, coherent_family_share=share,
                         coherent_family_size=int(fam.loc[f, "n_reactions"]), coherent_family_pct_orphan=float(fam.loc[f, "pct_orphan"]),
                         coherent_family_sample=fam.loc[f, "sample_reactions"][:200], candidate_score=share * float(fam.loc[f, "pct_orphan"]) / 100))
    return pd.DataFrame(rows), pd.DataFrame(avg)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--topk", required=True)
    ap.add_argument("--rxn-clusters", default=str(HERE / "rhea_reaction_clusters.csv"))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    fin, avg = summarize(pd.read_csv(a.topk), pd.read_csv(a.rxn_clusters))
    fin.to_csv(out / "protein_clusters_family_final.csv", index=False)
    avg.to_csv(out / "avg_top1_score.csv", index=False)
    print(f"{len(fin)} clusters -> {out}")


if __name__ == "__main__":
    main()
