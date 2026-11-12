#!/usr/bin/env python3
"""final_family_table.py reproduces the archived tables of the 33 dark clusters from the archived top-50 file.  python test_final_family_table.py"""
from pathlib import Path

import pandas as pd

import final_family_table as ft

HERE = Path(__file__).resolve().parent
import os
TOPK = Path(os.environ.get("HORIZYN_TOPK", HERE / "horizyn_full_topk.csv"))      # set HORIZYN_TOPK; the archived top-50 file (Zenodo: horizyn_dark_top50.csv; data/ after `make data`)


def test_synthetic():
    rxn = pd.DataFrame({"rhea_lr_id": [1, 2, 3, 4], "reaction_cluster": [7, 7, -1, 8], "has_enzyme": [False, True, False, True],
                        "is_sam_methylation": [False] * 4, "reaction_smiles": list("abcd")})
    topk = pd.DataFrame({"cluster": [1, 1, 1, 1], "accession": list("wxyz"), "rank": [1, 1, 1, 1], "rhea_lr_id": [1, 1, 3, 4], "score": [0.5, 0.7, 0.9, 0.1]})
    fin, avg = ft.summarize(topk, rxn)
    r = fin.iloc[0]
    assert r.coherent_dominant_family == 7 and r.coherent_family_share == 0.5 and r.pct_noise_votes == 25.0 and r.coherent_family_pct_orphan == 50.0
    assert abs(r.candidate_score - 0.25) < 1e-9 and avg.iloc[0].avg_top1_score == 0.55


def test_archived():
    topk = TOPK if TOPK.exists() else None
    if topk is None:
        print("archived top-50 file not present: skipped")
        return
    fin, avg = ft.summarize(pd.read_csv(topk), pd.read_csv(HERE / "rhea_reaction_clusters.csv"))
    old = pd.read_csv(HERE / "protein_clusters_family_final.csv")
    m = old.merge(fin, on="cluster", suffixes=("_old", "_new"))
    assert len(m) == len(old) == len(fin) == 33
    for col in ["n", "pct_noise_votes", "coherent_dominant_family", "coherent_family_share", "coherent_family_size", "coherent_family_pct_orphan", "coherent_family_sample"]:
        assert (m[f"{col}_old"] == m[f"{col}_new"]).all(), col
    assert ((m.candidate_score_old - m.candidate_score_new).abs() < 1e-5).all()
    oa = pd.read_csv(HERE / "avg_top1_score.csv").merge(avg, on="cluster", suffixes=("_old", "_new"))
    assert ((oa.avg_top1_score_old - oa.avg_top1_score_new).abs() < 1e-9).all()


if __name__ == "__main__":
    test_synthetic()
    test_archived()
    print("all tests passed")
