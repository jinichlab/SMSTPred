#!/usr/bin/env python3
"""
Stage 6 motif analysis -- step 2b: turn mmseqs clusters into per-class
representative sets.

Reads the mmseqs easy-cluster output of corpus_all.fasta (IDs are 'Entry|class').
Each cluster (a ~50%-identity family) gets one representative; the cluster's
class is the majority class of its members. Writes one FASTA of representatives
per class -- so downstream STREME sees one sequence per family, not thousands of
near-duplicate paralogs.

Outputs: cluster_reps.csv, reps_{Metabolite,Protein,RNA,DNA}.fasta

  python 2_make_reps.py [--dir <folder>]    the folder with clust_cluster.tsv / clust_rep_seq.fasta (default: this folder; a release: <release>/motif_analysis)
"""
import argparse
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

IO = Path(__file__).resolve().parent
LABELS = ["Metabolite", "Protein", "RNA", "DNA"]


def read_fasta(path):
    seqs, sid, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if sid is not None:
                    seqs[sid] = "".join(buf)
                sid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if sid is not None:
        seqs[sid] = "".join(buf)
    return seqs


def main():
    global IO
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=str(IO))
    IO = Path(ap.parse_args().dir)
    members = defaultdict(list)
    with open(IO / "clust_cluster.tsv") as fh:
        for line in fh:
            rep, mem = line.rstrip("\n").split("\t")
            members[rep].append(mem)
    rep_seq = read_fasta(IO / "clust_rep_seq.fasta")

    out = {l: [] for l in LABELS}
    rows = []
    for rep, mems in members.items():
        classes = [m.split("|")[1] for m in mems if "|" in m]
        if not classes:
            continue
        maj, n = Counter(classes).most_common(1)[0]
        purity = n / len(classes)
        entry = rep.split("|")[0]
        seq = rep_seq.get(rep)
        rows.append({"rep": entry, "class": maj, "cluster_size": len(mems),
                     "purity": round(purity, 3)})
        if maj in out and seq:
            out[maj].append((entry, seq))

    pd.DataFrame(rows).to_csv(IO / "cluster_reps.csv", index=False)
    print(f"{len(rows)} clusters (families)", flush=True)
    for lab in LABELS:
        with open(IO / f"reps_{lab}.fasta", "w") as fa:
            for e, s in out[lab]:
                fa.write(f">{e}\n{s}\n")
        print(f"  reps_{lab}.fasta: {len(out[lab])} family representatives", flush=True)


if __name__ == "__main__":
    main()
