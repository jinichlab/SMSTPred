#!/usr/bin/env python3
"""
Stage 1 -- MMseqs2 sequence clustering of the labelled sets, and the cluster-disjoint train / validation splits.

  python mmseqs_clusters.py sam --subset benchmark        # the SAM proteins with in_benchmark  (<= 1,000 aa): Fig 2A, Table 1, S2 / S3, CV
  python mmseqs_clusters.py sam --subset deployment       # the SAM proteins with in_deployment (<= 2,000 aa): deployment model, OOF thresholds
  python mmseqs_clusters.py up                            # UP enzymes (in_benchmark) + the SAM benchmark proteins: Fig 2B

Reads labels_SAM.csv / labels_UP.csv (substrate typing, 0_uniprot_pull/bin_substrates.py) for the Entry list, the four labels and the in_* flags, and the
sequences from sam_mtases.tsv.gz / up_enzymes.tsv.gz of the same release folder. No embedding is needed, so this step runs right after the typing and the
splits are fixed before anything is embedded. Replaces the ProtT5-cosine clustering of sam_clustering.py / up_clustering.py, which needed the embeddings first
and made one of the compared embedders define the split; those scripts stay for reproducing the archived paper.

MMseqs2 easy-cluster, the settings of the side project ~/Code/mmseqs_clustering that produced the MMseqs-split benchmark:
  --min-seq-id 0.3 (30 % identity) -c 0.8 --cov-mode 0 (80 % of BOTH sequences aligned) -s 7.5 (most sensitive) --cluster-mode 0 (greedy set cover) --cluster-reassign 1

Outputs (<release>/clusters/, tag = <set>_<subset>_id30):
  clusters_<tag>.csv          Entry, cluster, representative, cluster_size, the four labels (UP: also in_sam, sam_cluster)
  cluster_summary_<tag>.txt   cluster counts, singletons, the biggest clusters and their classes
  <tag>_cluster.tsv, <tag>_rep_seq.fasta, <tag>_params.txt, <tag>.fasta   the raw MMseqs2 output, the representatives, versions / parameters, the input
  splits_<tag>.json, splits_summary_<tag>.csv   (SAM only) 16 cluster-disjoint 90 / 10 train / validation splits, seed 42
and two columns written into the label table: cluster_<subset> and cluster_size_<subset> (empty for proteins outside the subset). Re-running
bin_substrates.py rewrites the label table without them: cluster after typing.

The split rule is the repo's split_clusters_train_val: shuffle the clusters, the first 90 % of CLUSTERS are train, the rest validation, redraw until every
label has a positive in validation; the seed-42 stream is consumed one iteration after another, so the first 6 splits are the same whether 6 or 16 are made.
Needs mmseqs on the PATH (or $MMSEQS); numpy and pandas.
"""
import argparse
import hashlib
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
PULL = Path(os.environ.get("UNIPROT_PULL_DIR", REPO / "data" / "uniprot_pull"))
LABELS = ["Metabolite", "Protein", "RNA", "DNA"]
FILES = {"sam": ("labels_SAM.csv", "sam_mtases.tsv.gz"), "up": ("labels_UP.csv", "up_enzymes.tsv.gz")}


def newest_release(pull_dir, release=None):
    if release:
        return Path(pull_dir) / f"release_{release}"
    rel = sorted(p for p in Path(pull_dir).glob("release_*") if p.is_dir() and not p.name.startswith("_"))
    if not rel:
        sys.exit(f"no release_* folder in {pull_dir}")
    return rel[-1]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def load_set(rel, which, flag):
    """Label rows of the proteins with the flag set, plus their sequences. Fails loudly on a missing or empty sequence."""
    lab_name, tsv_name = FILES[which]
    lab = pd.read_csv(rel / lab_name, low_memory=False)
    if lab["Entry"].duplicated().any():
        sys.exit(f"{lab_name}: Entry is not unique")
    if flag not in lab:
        sys.exit(f"{lab_name} has no column {flag}: run bin_substrates.py first")
    sub = lab[lab[flag].astype(bool)].reset_index(drop=True)
    seq = pd.read_csv(rel / tsv_name, sep="\t", usecols=["Entry", "Sequence"], low_memory=False).drop_duplicates("Entry").set_index("Entry")["Sequence"]
    sub["Sequence"] = sub["Entry"].map(seq).astype(str).str.strip().str.upper()
    bad = sub["Sequence"].isin(["", "NAN"])
    if bad.any():
        sys.exit(f"{int(bad.sum())} proteins have no sequence in {tsv_name}, e.g. {list(sub.loc[bad, 'Entry'][:5])}")
    return lab, sub


def write_fasta(path, entries, seqs):
    with open(path, "w") as f:
        for e, s in zip(entries, seqs):
            f.write(f">{e}\n{s}\n")


def run_mmseqs(fasta, prefix, p, threads, log):
    exe = os.environ.get("MMSEQS") or shutil.which("mmseqs")
    if not exe:
        sys.exit("mmseqs not found: put it on the PATH (bioconda: mmseqs2) or set $MMSEQS")
    version = subprocess.run([exe, "version"], capture_output=True, text=True).stdout.strip()
    tmp = tempfile.mkdtemp(prefix="mmseqs_", dir=os.environ.get("TMPDIR"))
    cmd = [exe, "easy-cluster", str(fasta), str(prefix), tmp, "--min-seq-id", str(p["min_seq_id"]), "-c", str(p["cov"]), "--cov-mode", str(p["cov_mode"]),
           "-s", str(p["sens"]), "--cluster-mode", str(p["cluster_mode"]), "--cluster-reassign", str(p["reassign"]), "--threads", str(threads)]
    try:
        with open(log, "w") as lf:
            subprocess.run(cmd, check=True, stdout=lf, stderr=subprocess.STDOUT)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return version


def read_clusters(tsv):
    clusters = {}
    for line in open(tsv):
        rep, mem = line.rstrip("\n").split("\t")
        clusters.setdefault(rep, []).append(mem)
    return clusters


def number_clusters(raw):
    """Stable ids: the biggest cluster is 0 (ties: representative Entry)."""
    order = sorted(raw.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    clusters = {cid: sorted(ms) for cid, (_, ms) in enumerate(order)}
    rep_of = {cid: rep for cid, (rep, _) in enumerate(order)}
    return clusters, rep_of, {e: cid for cid, ms in clusters.items() for e in ms}


def split_clusters_train_val(clusters, y, train_ratio=0.9, max_attempts=1000):
    """The repo's rule (2_training/mlp_embedding_comparison.py): shuffle the cluster ids, first 90 % train, redraw until every label is in validation."""
    cids = list(clusters.keys())
    for _ in range(max_attempts):
        random.shuffle(cids)
        k = int(train_ratio * len(cids))
        tr = [i for c in cids[:k] for i in clusters[c]]
        va = [i for c in cids[k:] for i in clusters[c]]
        if all(y[va][:, j].sum() > 0 for j in range(y.shape[1])):
            return tr, va, cids[k:]
    raise ValueError("no valid cluster split")


def make_splits(entries, y, clusters, n_iter, seed):
    idx = {e: i for i, e in enumerate(entries)}
    cl_idx = {cid: [idx[e] for e in ms] for cid, ms in clusters.items()}
    random.seed(seed)
    iterations, rows = [], []
    for it in range(n_iter):
        tr, va, val_cids = split_clusters_train_val(cl_idx, y, 0.9)
        trs, vas = {entries[i] for i in tr}, {entries[i] for i in va}
        assert not trs & vas
        assert all(len({e in vas for e in ms}) == 1 for ms in clusters.values()), "a cluster straddles the split"
        iterations.append({"train": sorted(trs), "val": sorted(vas)})
        pos = y[va].sum(0).astype(int)
        rows.append({"iteration": it + 1, "n_train": len(tr), "n_val": len(va), "val_fraction": round(len(va) / len(entries), 4), "n_val_clusters": len(val_cids),
                     "largest_val_cluster": max(len(clusters[c]) for c in val_cids), **{f"val_{l}": int(p) for l, p in zip(LABELS, pos)}})
    return iterations, pd.DataFrame(rows)


def summary(tag, entries, clusters, rep_of, lab_by_entry):
    sizes = np.array(sorted((len(m) for m in clusters.values()), reverse=True))
    lines = [f"tag: {tag}", f"proteins: {len(entries)}", f"clusters: {len(clusters)}",
             f"singletons: {int((sizes == 1).sum())} ({(sizes == 1).mean():.1%} of clusters, {(sizes == 1).sum() / len(entries):.1%} of proteins)",
             f"size median / mean / max: {int(np.median(sizes))} / {sizes.mean():.2f} / {sizes[0]}", f"ten largest clusters: {sizes[:10].tolist()}",
             f"proteins in the 10 largest clusters: {int(sizes[:10].sum())} ({sizes[:10].sum() / len(entries):.1%})"]
    for cid in range(min(5, len(clusters))):
        lines.append(f"  cluster {cid} (rep {rep_of[cid]}): size {len(clusters[cid])}, label counts {lab_by_entry.loc[clusters[cid], LABELS].sum().astype(int).to_dict()}")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("which", choices=["sam", "up"])
    ap.add_argument("--subset", choices=["benchmark", "deployment"], default="benchmark", help="sam only; the UP clustering uses the benchmark set")
    ap.add_argument("--pull-dir", default=str(PULL))
    ap.add_argument("--release")
    ap.add_argument("--out", help="output folder (default <release>/clusters)")
    ap.add_argument("--min-seq-id", type=float, default=0.3)
    ap.add_argument("--cov", type=float, default=0.8)
    ap.add_argument("--cov-mode", type=int, default=0)
    ap.add_argument("--sens", type=float, default=7.5)
    ap.add_argument("--cluster-mode", type=int, default=0)
    ap.add_argument("--reassign", type=int, default=1)
    ap.add_argument("--tag-suffix", default="")
    ap.add_argument("--n-iterations", type=int, default=16)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    a = ap.parse_args()
    if a.which == "up" and a.subset != "benchmark":
        sys.exit("the UP clustering is made for the benchmark set only (Fig 2B)")

    rel = newest_release(a.pull_dir, a.release)
    out = Path(a.out) if a.out else rel / "clusters"
    out.mkdir(parents=True, exist_ok=True)
    flag = f"in_{a.subset}"
    tag = f"{a.which}_{a.subset}_id{int(round(a.min_seq_id * 100))}{a.tag_suffix}"
    print(f"{tag}: release {rel.name}, flag {flag}, {a.threads} threads", flush=True)

    full_lab, sub = load_set(rel, a.which, flag)
    sam_only = pd.DataFrame()
    if a.which == "up":  # the SAM benchmark proteins go into the same clustering, so UP homologues of a held-out SAM protein can be found
        _, sam = load_set(rel, "sam", "in_benchmark")
        sam_only = sam[~sam["Entry"].isin(sub["Entry"])]
        sub = pd.concat([sub.assign(in_sam=sub["Entry"].isin(sam["Entry"]).astype(int)), sam_only.assign(in_sam=1)], ignore_index=True)
        print(f"UP benchmark proteins {len(sub) - len(sam_only)} + SAM benchmark proteins not in UP {len(sam_only)}", flush=True)
    entries = list(sub["Entry"])
    fasta = out / f"{tag}.fasta"
    write_fasta(fasta, entries, sub["Sequence"])
    lens = sub["Sequence"].str.len()
    print(f"{len(entries)} sequences, length min/median/max = {lens.min()}/{int(lens.median())}/{lens.max()}", flush=True)

    params = {"min_seq_id": a.min_seq_id, "cov": a.cov, "cov_mode": a.cov_mode, "sens": a.sens, "cluster_mode": a.cluster_mode, "reassign": a.reassign}
    prefix = out / tag
    version = run_mmseqs(fasta, prefix, params, a.threads, out / f"{tag}_mmseqs.log")
    raw = read_clusters(f"{prefix}_cluster.tsv")
    members = [m for ms in raw.values() for m in ms]
    assert sorted(members) == sorted(entries), f"the cluster file covers {len(set(members))} proteins, the input has {len(entries)}"
    clusters, rep_of, cid_of = number_clusters(raw)

    table = pd.DataFrame({"Entry": entries, "cluster": [cid_of[e] for e in entries]})
    table["representative"] = table["cluster"].map(rep_of)
    table["cluster_size"] = table["cluster"].map({c: len(m) for c, m in clusters.items()})
    if a.which == "up":
        table["in_sam"] = sub["in_sam"].values
        sam_clusters = out / "clusters_sam_benchmark_id30.csv"
        if sam_clusters.exists():
            table["sam_cluster"] = table["Entry"].map(pd.read_csv(sam_clusters).set_index("Entry")["cluster"])
        else:
            print("  (no SAM-only benchmark clustering yet, so sam_cluster stays empty: run `sam --subset benchmark` first)")
    table = pd.concat([table, sub[LABELS].reset_index(drop=True)], axis=1).sort_values(["cluster", "Entry"])
    table.to_csv(out / f"clusters_{tag}.csv", index=False)
    text = summary(tag, entries, clusters, rep_of, sub.set_index("Entry"))
    print(text)
    (out / f"cluster_summary_{tag}.txt").write_text(text + "\n")
    (out / f"{tag}_params.txt").write_text(f"mmseqs {version}\ndate {datetime.now().astimezone().isoformat(timespec='seconds')}\nrelease {rel.name}\nflag {flag}\n"
                                             f"parameters {params}\ninput {fasta.name} ({len(entries)} sequences, sha256 {sha256(fasta)})\n")

    if a.which == "sam":
        y = sub[LABELS].values.astype(np.float32)
        iterations, summ = make_splits(entries, y, clusters, a.n_iterations, a.seed)
        (out / f"splits_{tag}.json").write_text(json.dumps({"tag": tag, "seed": a.seed, "n_proteins": len(entries), "n_clusters": len(clusters), "iterations": iterations}))
        summ.to_csv(out / f"splits_summary_{tag}.csv", index=False)
        print(f"\nwrote splits_{tag}.json\n{summ.to_string(index=False)}\nlabel totals in the set: {dict(zip(LABELS, y.sum(0).astype(int)))}")

    # the cluster columns go into the label table (atomically; proteins outside the subset stay empty)
    lab_name = FILES[a.which][0]
    cols = [f"cluster_{a.subset}", f"cluster_size_{a.subset}"]
    full_lab = full_lab.drop(columns=[c for c in cols if c in full_lab])
    by = table.drop_duplicates("Entry").set_index("Entry")
    full_lab[cols[0]] = full_lab["Entry"].map(by["cluster"]).astype("Int64")
    full_lab[cols[1]] = full_lab["Entry"].map(by["cluster_size"]).astype("Int64")
    tmp = rel / f".{lab_name}.tmp"
    full_lab.to_csv(tmp, index=False)
    os.replace(tmp, rel / lab_name)
    print(f"\nadded {cols} to {lab_name}; outputs in {out}")


if __name__ == "__main__":
    main()
