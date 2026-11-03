#!/usr/bin/env python3
"""The per-organism table behind Fig 5 and 4_phylo_trends/lineage_statistics.py, built from a refresh release.

  python build_proteome_vs_sam_table.py [--release <folder>]      (default: $RELEASE_DIR, else the newest release_* of $UNIPROT_PULL_DIR)

One row per selected organism that has a proteome in the release (data/selected_1259_organisms.csv joined by Proteome Id to <release>/proteome_summary.csv):
  organism       the proteome's file stem (as in the archived table)
  domain         Bacteria / Eukaryota / Archaea, from the taxonomic lineage
  proteome_size  proteins in the reference proteome (n_all)
  n_sam          proteins with the SAM-dependent MTase superfamily domain IPR029063 (every length)
  n_predicted    of those, the ones <= 2,000 aa that the ensemble scored
  n_DNA, n_RNA, n_Protein, n_Metabolite   the ensemble's calls among the scored proteins (ALL calls, in domain or not; a protein with two calls counts in both)
so a substrate count misses the few IPR029063 proteins longer than 2,000 aa (n_sam - n_predicted). Writes <release>/fig5/proteome_vs_sam_bylabel.csv.
"""
import argparse
import os
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
META = Path(os.environ.get("SMSTPRED_DATA", REPO / "data")) / "selected_1259_organisms.csv"


def domain_of(lineage):
    lineage = lineage if isinstance(lineage, str) else ""
    return "Bacteria" if "Bacteria" in lineage else "Eukaryota" if "Eukaryota" in lineage else "Archaea" if "Archaea" in lineage else "other"


def build(rel, meta):
    summ = pd.read_csv(rel / "proteome_summary.csv")
    summ = summ[summ.status == "ok"][["proteome_id", "stem", "n_all", "n_ipr_all"]]
    pred = pd.read_csv(rel / "predictions_all_proteomes.csv.gz", usecols=["proteome", "call_DNA", "call_RNA", "call_Protein", "call_Metabolite"])
    g = pred.groupby("proteome").agg(n_predicted=("call_DNA", "size"), n_DNA=("call_DNA", "sum"), n_RNA=("call_RNA", "sum"),
                                     n_Protein=("call_Protein", "sum"), n_Metabolite=("call_Metabolite", "sum"))
    d = meta.merge(summ, left_on="Proteome Id", right_on="proteome_id", how="inner").merge(g, left_on="stem", right_index=True, how="left")
    d[["n_predicted", "n_DNA", "n_RNA", "n_Protein", "n_Metabolite"]] = d[["n_predicted", "n_DNA", "n_RNA", "n_Protein", "n_Metabolite"]].fillna(0).astype(int)
    out = pd.DataFrame({"organism": d["stem"], "domain": d["Taxonomic lineage"].map(domain_of), "proteome_size": d["n_all"], "n_sam": d["n_ipr_all"],
                        "n_predicted": d["n_predicted"], "n_DNA": d["n_DNA"], "n_RNA": d["n_RNA"], "n_Protein": d["n_Protein"], "n_Metabolite": d["n_Metabolite"]})
    return out.sort_values("organism").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release")
    a = ap.parse_args()
    sys.path.insert(0, str(REPO / "1_clustering"))
    from mmseqs_clusters import newest_release
    rel = Path(a.release or os.environ.get("RELEASE_DIR") or newest_release(Path(os.environ.get("UNIPROT_PULL_DIR", REPO / "data" / "uniprot_pull")), os.environ.get("RELEASE")))
    out = build(rel, pd.read_csv(META))
    (rel / "fig5").mkdir(exist_ok=True)
    out.to_csv(rel / "fig5" / "proteome_vs_sam_bylabel.csv", index=False)
    print(f"{len(out)} organisms ({out.domain.value_counts().to_dict()}) -> {rel}/fig5/proteome_vs_sam_bylabel.csv")


if __name__ == "__main__":
    main()
