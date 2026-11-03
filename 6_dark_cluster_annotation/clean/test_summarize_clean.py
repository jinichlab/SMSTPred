#!/usr/bin/env python3
"""summarize_clean_clusters.py reproduces the tier, distance, top EC and MTase call of the 33 archived dark clusters.  python test_summarize_clean.py"""
from pathlib import Path

import pandas as pd

import summarize_clean_clusters as sc

HERE = Path(__file__).resolve().parent
PVALUE = Path("/home/nroethler/Code/abhiram/CLEAN/app/results/inputs/dark33_pvalue.csv")      # CLEAN's raw p-value output of the archived run (Zenodo / the CLEAN checkout)


def test_tier():
    assert [sc.tier(x) for x in (5.99, 6.0, 7.99, 8.0)] == ["reliable", "marginal", "marginal", "unreliable"]


def test_synthetic():
    maxsep = {"a": [("2.1.1.1", .9)], "b": [("2.1.1.1", .8)], "c": [("2.5.1.16", .7)]}
    pvalue = {"a": [("2.1.1.1", 5.0)], "b": [("2.1.1.1", 7.0)], "c": [("2.5.1.16", 9.0)]}
    m = pd.DataFrame({"Entry": ["a", "b", "c"], "cluster": [1, 1, 2], "lineage": ["Plant"] * 3})
    t = sc.summarize(maxsep, pvalue, m, {"2.1.1.1": "nicotinamide N-methyltransferase"}).set_index("cluster")
    assert t.loc[1, "top_ec"] == "EC:2.1.1.1" and t.loc[1, "d_median"] == 6.0 and t.loc[1, "confidence"].startswith("marginal") and t.loc[1, "mtase"] == "yes"
    assert t.loc[2, "mtase"] == "No" and t.loc[2, "confidence"].startswith("unreliable") and t.loc[1, "enzyme"] == "nicotinamide N-methyltransferase"


def test_archived():
    if not PVALUE.exists():
        print("archived CLEAN p-value output not present: skipped")
        return
    t = sc.summarize(sc.read_calls(HERE / "dark33_maxsep.csv"), sc.read_calls(PVALUE), pd.read_csv(HERE.parent / "dark33_members.csv")).set_index("cluster")
    old = pd.read_csv(HERE / "clean_dark33_cluster_summary.csv").set_index("cluster")
    per = pd.read_csv(HERE / "clean_dark33_per_cluster_ec.csv").set_index("cluster")
    assert len(t) == len(old) == 33
    old_ec = per.top_ec.reindex(t.index)
    differ = t.index[t.top_ec != old_ec].tolist()
    assert differ in ([], [733]), differ                      # 733: the archived parser cut the preliminary number EC:2.1.1.n7 to "EC:2.1.1." (fixed in parse_clean_results.py)
    if differ:
        assert t.loc[733, "top_ec"] == "EC:2.1.1.n7" and old_ec.loc[733] == "EC:2.1.1."
    assert (t.confidence.str.split().str[0] == old.confidence.str.split().str[0].reindex(t.index)).all(), "tier"
    d_old = old.confidence.str.extract(r"d=([\d.]+)")[0].astype(float).reindex(t.index)
    assert ((t.d_median - d_old).abs() < 0.06).all(), "distance"
    assert (t.mtase == old.mtase.reindex(t.index).str.strip().replace({"yes": "yes", "No": "No", "no": "No"})).all(), "mtase"


if __name__ == "__main__":
    test_tier()
    test_synthetic()
    test_archived()
    print("all tests passed")
