#!/usr/bin/env python3
"""Tests of cluster_selection.py: the rules on synthetic clusters, and (when the archived files are present) that they recover the archived Fig 7 / Table S1-S2 clusters.  python test_cluster_selection.py"""
import re
from pathlib import Path

import pandas as pd

import cluster_selection as cs

HERE = Path(__file__).resolve().parent
REPO = HERE.parent


def synthetic():
    rows = []

    def add(cluster, n, substrate, org, reviewed=0):
        for i in range(n):
            rows.append({"Entry": f"c{cluster}_{i}", "Organism": org, "Substrate": substrate, "cluster": cluster, "Reviewed": "reviewed" if i < reviewed else "unreviewed", "Protein names": "x"})
    add(1, 20, "Metabolite", "act", 0)                                    # dark Actinomycetota
    add(2, 20, "Metabolite", "act", 1)                                    # same but one reviewed member -> not dark
    add(3, 18, "Metabolite", "act"); add(3, 2, "RNA", "act")              # 90 % Actinomycetota and Metabolite -> dark
    add(4, 16, "Metabolite", "act"); add(4, 4, "Metabolite", "other")     # 80 % Actinomycetota & Metabolite exactly -> dark
    add(5, 15, "Metabolite", "act"); add(5, 5, "RNA", "act")              # lineage 100 % but Metabolite 75 % -> not dark
    add(6, 20, "Protein", "met", 2)                                       # metazoan protein cluster
    add(7, 20, "Protein", "met", 0)                                       # no anchor
    add(8, 40, "DNA", "vert", 1)                                          # vertebrate DNA cluster
    add(9, 20, "DNA", "vert", 1)                                          # too few vertebrates
    add(10, 22, "DNA", "vert", 0); add(10, 10, "RNA", "met", 0)          # 22 vertebrates (< 30): selected only with a reviewed vertebrate DNMT in it
    meta = pd.DataFrame({"Organism": ["act", "other", "met", "vert"], "Taxonomic lineage": ["Bacteria, Actinomycetota", "Bacteria, Bacillota", "Eukaryota, Metazoa", "Eukaryota, Metazoa, Vertebrata"],
                         "is_metazoa": [False, False, True, True], "is_vertebrate": [False, False, False, True]})
    return cs.annotate(pd.DataFrame(rows).drop(columns=["Reviewed", "Protein names"]), pd.DataFrame(rows)[["Entry", "Reviewed", "Protein names"]], meta)


def test_rules():
    t = cs.cluster_table(synthetic())
    assert set(cs.select_dark(t).cluster) == {1, 3, 4}
    assert set(cs.select_metazoa_protein(t).cluster) == {6}
    assert set(cs.select_vertebrate_dna(t).cluster) == {8}, "no anchor, so cluster 10 is not selected"
    d = synthetic()
    first = (d.cluster == 10) & (d.Entry == "c10_0")
    d.loc[first, "reviewed"] = True
    d.loc[first, "Protein names"] = "DNA (cytosine-5)-methyltransferase 3A"
    d["rev_vert_dnmt"] = d.reviewed & d.vertebrate & d["Protein names"].fillna("").str.lower().str.contains("dna")
    assert set(cs.select_vertebrate_dna(cs.cluster_table(d)).cluster) == {8, 10}, "the anchor clause adds a mostly-DNA vertebrate cluster with a reviewed DNMT"


def test_crosswalk():
    new = pd.DataFrame({"Entry": list("abcdefgh"), "cluster": [0, 0, 0, 0, 1, 1, 1, -1]})
    old = pd.DataFrame({"Entry": list("abcdefg"), "cluster": [5, 5, 5, 6, 6, 6, 6]})
    cw = cs.crosswalk(new, old, {"s": [5, 6]}).set_index("old_cluster")
    assert cw.loc[5, "new_cluster"] == 0 and cw.loc[5, "shared"] == 3 and abs(cw.loc[5, "jaccard"] - 3 / 4) < 1e-9
    assert cw.loc[6, "new_cluster"] == 1 and cw.loc[6, "shared"] == 3 and abs(cw.loc[6, "old_in_new"] - 0.75) < 1e-9


def test_recovers_archived():
    coords, rev = HERE / "sam_mtase_umap_prostt5_coords.csv", HERE / "reviewed_status_full.csv"
    if not (coords.exists() and rev.exists() and cs.META_CSV.exists()):
        print("archived files not present: skipped")
        return
    d = cs.annotate(pd.read_csv(coords, usecols=["Entry", "Organism", "Substrate", "cluster"]), pd.read_csv(rev, usecols=["Entry", "Reviewed", "Protein names"]), pd.read_csv(cs.META_CSV))
    t = cs.cluster_table(d)
    dark = cs.select_dark(t)
    assert len(dark) == 33 and (dark.lineage == "Actinomycetota").sum() == 14 and (dark.lineage == "Plant").sum() == 19, "the 33 dark clusters of the archived map"
    assert set(dark.cluster) == set(pd.read_csv(REPO / "6_dark_cluster_annotation" / "dark33_members.csv").cluster)
    assert set(cs.select_metazoa_protein(t).cluster) == set(cs.ARCHIVED_A), "Fig 7A"
    assert set(cs.select_vertebrate_dna(t).cluster) == set(cs.ARCHIVED_B), "Fig 7B"


if __name__ == "__main__":
    test_rules()
    test_crosswalk()
    test_recovers_archived()
    print("all cluster-selection tests passed")
