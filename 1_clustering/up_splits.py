#!/usr/bin/env python3
"""
The UP training / validation sets of the SAM-vs-UP comparison, one pair per SAM test set of the bake-off.

  python up_splits.py [--release <release folder>] [--val-frac 0.15] [--seed 42]     (defaults: ../refresh_config.json; the release: $RELEASE_DIR, else the newest release_*)

The question: does a model trained on the broad UP enzyme set predict the substrate type of held-out SAM methyltransferases as well as one trained on SAM? Both are scored on the SAME test
proteins: the 20 % held-out SAM clusters of split i (tvt_splits.py). For the UP model, for every split i:
  1. every UP cluster that contains a test protein of split i is REMOVED from the UP set (mmseqs_clusters.py up clusters UP and SAM together at 30 % identity, so this takes out the test proteins
     themselves and every UP protein that is a homologue of one: the UP model gets no homology advantage over the cluster-disjoint SAM model);
  2. of the remaining UP clusters, 15 % (round(val_frac * n), at random) are the VALIDATION set (early stopping and the best weights, as for SAM), the rest the TRAIN set.
The UP set contains the SAM proteins too (the 7,553 benchmark proteins are in labels_UP.csv): those of SAM's own train / validation clusters stay in the UP pool, those whose cluster was
removed go with it. The random draw of split i is seeded by (seed, i), so it does not depend on how many splits are made.

Reads <release>/clusters/clusters_up_benchmark_id30.csv (Entry, cluster, the four labels) and <release>/clusters/splits_tvt_sam_benchmark_id30_t20_v15.json; writes
  <release>/clusters/splits_tvt_up_id30_t20_v15.json   {"iterations": [{"train": [...], "validation": [...], "removed_clusters": [...], "n_removed_proteins": n}, ...]}
  <release>/clusters/splits_tvt_summary_up_id30_t20_v15.csv   per split and partition (train / validation / removed): clusters, proteins, positives of the four labels
(t20 = the SAM test fraction the splits were made for, v15 = the UP validation fraction.)
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from mmseqs_clusters import newest_release

LABELS = ["Metabolite", "Protein", "RNA", "DNA"]
REPO = Path(__file__).resolve().parents[1]
CFG = json.loads((REPO / "refresh_config.json").read_text()) if (REPO / "refresh_config.json").exists() else {}


def make_up_split(up, test_entries, val_frac, seed, split_index):
    """up: table with Entry and cluster. -> (train clusters, validation clusters, removed clusters) as sets of cluster ids."""
    removed = set(up.loc[up.Entry.isin(set(test_entries)), "cluster"])
    remaining = np.array(sorted(set(up.cluster) - removed))
    rng = np.random.default_rng(np.random.SeedSequence([seed, 2, split_index]))
    perm = rng.permutation(remaining)
    n_val = int(round(val_frac * len(remaining)))
    assert 0 < n_val < len(remaining)
    return set(perm[n_val:].tolist()), set(perm[:n_val].tolist()), removed


def summarize(up, parts_by_split):
    rows = []
    for i, (tr, va, rm) in enumerate(parts_by_split, 1):
        for name, cl in (("train", tr), ("validation", va), ("removed", rm)):
            sub = up[up.cluster.isin(cl)]
            rows.append({"split": i, "partition": name, "clusters": len(cl), "proteins": len(sub), **{l: int(sub[l].sum()) for l in LABELS}})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release")
    ap.add_argument("--val-frac", type=float, default=CFG.get("splits", {}).get("val_frac", 0.15))
    ap.add_argument("--test-frac", type=float, default=CFG.get("splits", {}).get("test_frac", 0.2), help="only to find the SAM split file")
    ap.add_argument("--seed", type=int, default=CFG.get("seed", 42))
    a = ap.parse_args()
    rel = Path(a.release or os.environ.get("RELEASE_DIR") or newest_release(Path(os.environ.get("UNIPROT_PULL_DIR", REPO / "data" / "uniprot_pull")), os.environ.get("RELEASE")))
    tag = f"t{round(a.test_frac * 100)}_v{round(a.val_frac * 100)}"
    sam_file = rel / "clusters" / f"splits_tvt_sam_benchmark_id30_t{round(a.test_frac * 100)}_v{round(CFG.get('splits', {}).get('val_frac', 0.15) * 100)}.json"
    up_file = rel / "clusters" / "clusters_up_benchmark_id30.csv"
    for f in (sam_file, up_file):
        if not f.exists():
            sys.exit(f"{f} not found: run 1_clustering/run_mmseqs_clusters.sh (CLUSTER_WHAT=up for the UP clustering; the SAM benchmark clustering writes the splits)")
    up = pd.read_csv(up_file, usecols=["Entry", "cluster", "in_sam", *LABELS], low_memory=False)
    sam_iters = json.loads(sam_file.read_text())["iterations"]
    ent_by_cluster = up.groupby("cluster").Entry.apply(list).to_dict()
    pack = lambda cl: sorted(e for c in cl for e in ent_by_cluster[c])
    parts, out_iters = [], []
    for i, it in enumerate(sam_iters, 1):
        tr, va, rm = make_up_split(up, it["test"], a.val_frac, a.seed, i)
        assert not (tr & va) and not (tr & rm) and not (va & rm) and len(tr | va | rm) == up.cluster.nunique()
        test = set(it["test"])
        assert not (test & set(pack(tr))) and not (test & set(pack(va))), "a SAM test protein is in the UP training or validation set"
        parts.append((tr, va, rm))
        out_iters.append({"train": pack(tr), "validation": pack(va), "removed_clusters": sorted(int(c) for c in rm), "n_removed_proteins": int(sum(len(ent_by_cluster[c]) for c in rm))})
    (rel / "clusters" / f"splits_tvt_up_id30_{tag}.json").write_text(json.dumps({
        "tag": f"up_id30_{tag}", "seed": a.seed, "val_frac": a.val_frac, "n_proteins": int(len(up)), "n_clusters": int(up.cluster.nunique()),
        "rule": "UP clusters containing a SAM test protein are removed; of the remaining UP clusters val_frac are validation, the rest train", "iterations": out_iters}))
    summ = summarize(up, parts)
    summ.to_csv(rel / "clusters" / f"splits_tvt_summary_up_id30_{tag}.csv", index=False)

    print(f"{len(up)} UP proteins in {up.cluster.nunique()} clusters; test sets from {sam_file.name}\n")
    fmt = lambda r: f"{r.proteins:6d} in {r.clusters:5d} cl ({r.Metabolite:6d}/{r.Protein:5d}/{r.RNA:5d}/{r.DNA:4d})"
    print(f"{'split':>5s} | {'TRAIN  proteins in clusters (Met/Prot/RNA/DNA)':^50s} | {'VALIDATION':^50s} | {'REMOVED (clusters with a test protein)':^45s}")
    for i in range(1, len(parts) + 1):
        s = summ[summ.split == i].set_index("partition")
        print(f"{i:5d} | {fmt(s.loc['train']):^50s} | {fmt(s.loc['validation']):^50s} | {fmt(s.loc['removed']):^45s}")
    v = summ[summ.partition == "validation"]
    print(f"\nfewest validation positives in any split: {({l: int(v[l].min()) for l in LABELS})}; validation proteins {int(v.proteins.min())}-{int(v.proteins.max())}")
    ins = up.set_index("Entry").in_sam.astype(bool)
    up_only = [int((~ins.loc[[e for c in rm for e in ent_by_cluster[c]]]).sum()) for _, _, rm in parts]
    sam_other = [int(out_iters[i]["n_removed_proteins"] - len(sam_iters[i]["test"]) - up_only[i]) for i in range(len(parts))]
    print(f"removed with the test clusters, per split: UP-only proteins (homologues of a test protein) {up_only}\n"
          f"                                          SAM train / validation proteins that share a joint cluster with a test protein {sam_other}")
    print(f"\nwrote {rel / 'clusters'}/splits_tvt_up_id30_{tag}.json and splits_tvt_summary_up_id30_{tag}.csv")


if __name__ == "__main__":
    main()
