#!/usr/bin/env python3
"""
Cluster-disjoint TRAIN / VALIDATION / TEST splits of the SAM benchmark proteins: 70 % / 10 % / 20 % of the CLUSTERS, drawn at random, 16 times.

  python tvt_splits.py [--release <release folder>] [--n 16] [--test-frac 0.2] [--val-frac 0.15] [--seed 42]     (defaults: ../refresh_config.json; the release: $RELEASE_DIR, else the newest release_* of $UNIPROT_PULL_DIR)

Reads <release>/labels_SAM.csv (the in_benchmark proteins with their cluster_benchmark, written by mmseqs_clusters.py sam --subset benchmark) and writes
  <release>/clusters/splits_tvt_sam_benchmark_id30_t20_v15.json   {"iterations": [{"train": [...], "validation": [...], "test": [...]}, ...]}  (UniProt entries)
  <release>/clusters/splits_tvt_summary_sam_benchmark_id30_t20_v15.csv   one row per split and partition: clusters, proteins and the four label counts
(t20_v15 = 20 % of the clusters test, 15 % validation; the fractions are part of the file names so that different choices never overwrite each other)
and prints the same table. The keys are "validation" and "test" on purpose: the 90 / 10 files of mmseqs_clusters.py call their one held-out part "val".

The rule: for every split the cluster ids are shuffled (numpy default_rng(seed), the stream is consumed one split after another, so the first k splits do not depend on
how many are made); the first round(test_frac * n_clusters) clusters are the test set, the next round(val_frac * n_clusters) the validation set, the rest train. Nothing else:
there is NO redraw and NO balancing by label or protein count, so what a partition contains is what chance gives (cluster sizes are very uneven: 248 of the 502 benchmark
clusters are singletons, the largest holds 830 proteins, so 20 % of the clusters is not 20 % of the proteins). The table is there to show whether every partition has enough
positives of every label; --min-pos N lists the splits where one has fewer than N.
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
TAG = "sam_benchmark_id30"
REPO = Path(__file__).resolve().parents[1]
CFG = json.loads((REPO / "refresh_config.json").read_text()) if (REPO / "refresh_config.json").exists() else {}   # the pipeline's one config file


def make_tvt_splits(clusters, n_splits=16, test_frac=0.2, val_frac=0.1, seed=42):
    """clusters: sorted array of the distinct cluster ids. -> [(train_clusters, val_clusters, test_clusters), ...] (sets of ids)."""
    clusters = np.asarray(sorted(set(clusters)))
    n = len(clusters)
    n_test, n_val = int(round(test_frac * n)), int(round(val_frac * n))
    assert n_test > 0 and n_val > 0 and n_test + n_val < n, (n, n_test, n_val)
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n_splits):
        perm = rng.permutation(clusters)
        out.append((set(perm[n_test + n_val:].tolist()), set(perm[n_test:n_test + n_val].tolist()), set(perm[:n_test].tolist())))
    return out


def summarize(df, splits):
    """One row per split x partition: clusters, proteins and positives of every label (df: Entry, cluster, the four labels)."""
    rows = []
    for i, parts in enumerate(splits, 1):
        for name, cl in zip(("train", "validation", "test"), parts):
            sub = df[df.cluster.isin(cl)]
            rows.append({"split": i, "partition": name, "clusters": len(cl), "proteins": len(sub), **{l: int(sub[l].sum()) for l in LABELS}})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", help="release folder (default: the newest release_* in $UNIPROT_PULL_DIR, default <repo>/data/uniprot_pull)")
    ap.add_argument("--n", type=int, default=CFG.get("splits", {}).get("n_splits", 16))
    ap.add_argument("--test-frac", type=float, default=CFG.get("splits", {}).get("test_frac", 0.2))
    ap.add_argument("--val-frac", type=float, default=CFG.get("splits", {}).get("val_frac", 0.15))
    ap.add_argument("--seed", type=int, default=CFG.get("seed", 42))
    ap.add_argument("--min-pos", type=int, default=10, help="list the splits where a validation or test partition has fewer than this many positives of a label")
    a = ap.parse_args()

    rel = Path(a.release or os.environ.get("RELEASE_DIR") or newest_release(Path(os.environ.get("UNIPROT_PULL_DIR", REPO / "data" / "uniprot_pull")), os.environ.get("RELEASE")))
    d = pd.read_csv(rel / "labels_SAM.csv", usecols=["Entry", *LABELS, "in_benchmark", "cluster_benchmark"])
    d = d[d.in_benchmark == True].rename(columns={"cluster_benchmark": "cluster"})
    if d.cluster.isna().any():
        sys.exit("some benchmark proteins have no cluster: run mmseqs_clusters.py sam --subset benchmark first")
    d["cluster"] = d.cluster.astype(int)
    splits = make_tvt_splits(d.cluster.unique(), a.n, a.test_frac, a.val_frac, a.seed)
    for tr, va, te in splits:
        assert not (tr & va) and not (tr & te) and not (va & te) and len(tr | va | te) == d.cluster.nunique()      # every cluster in exactly one partition
    out = rel / "clusters"
    out.mkdir(exist_ok=True)
    tag = f"{TAG}_t{round(a.test_frac * 100)}_v{round(a.val_frac * 100)}"
    ent = d.groupby("cluster").Entry.apply(list).to_dict()
    pack = lambda cl: sorted(e for c in cl for e in ent[c])
    (out / f"splits_tvt_{tag}.json").write_text(json.dumps({
        "tag": TAG, "seed": a.seed, "n_proteins": int(len(d)), "n_clusters": int(d.cluster.nunique()), "test_frac": a.test_frac, "val_frac": a.val_frac,
        "rule": "random clusters: first round(test_frac*n) = test, next round(val_frac*n) = validation, rest train; no redraw, no balancing",
        "iterations": [{"train": pack(tr), "validation": pack(va), "test": pack(te)} for tr, va, te in splits]}))
    summ = summarize(d, splits)
    summ.to_csv(out / f"splits_tvt_summary_{tag}.csv", index=False)

    print(f"{len(d)} proteins in {d.cluster.nunique()} clusters; per split {len(splits[0][2])} test, {len(splits[0][1])} validation, {len(splits[0][0])} train clusters (seed {a.seed})\n")
    fmt = lambda r: f"{r.proteins:5d} ({r.Metabolite:4d}/{r.Protein:4d}/{r.RNA:4d}/{r.DNA:3d})"
    print(f"{'split':>5s} | {'TRAIN  n (Met/Prot/RNA/DNA)':^32s} | {'VALIDATION':^32s} | {'TEST':^32s}")
    for i in range(1, len(splits) + 1):
        s = summ[summ.split == i].set_index("partition")
        print(f"{i:5d} | {fmt(s.loc['train']):^32s} | {fmt(s.loc['validation']):^32s} | {fmt(s.loc['test']):^32s}")
    v, t = summ[summ.partition == "validation"], summ[summ.partition == "test"]
    print("\nfewest positives in any split, validation / test:", {l: f"{int(v[l].min())} / {int(t[l].min())}" for l in LABELS})
    print("fewest proteins in any split, validation / test:", f"{int(v.proteins.min())} / {int(t.proteins.min())}", "| most:", f"{int(v.proteins.max())} / {int(t.proteins.max())}")
    low = summ[summ.partition.isin(["validation", "test"])]
    flag = [(int(r.split), r.partition, l, int(r[l])) for _, r in low.iterrows() for l in LABELS if r[l] < a.min_pos]
    print(f"\nvalidation/test partitions with fewer than {a.min_pos} positives of a label: {len(flag)} of {len(low) * len(LABELS)} (split, partition, label)")
    for f in flag:
        print("   split %d %-10s %-10s %d" % f)
    print(f"\nwrote {out}/splits_tvt_{tag}.json and splits_tvt_summary_{tag}.csv")


if __name__ == "__main__":
    main()
