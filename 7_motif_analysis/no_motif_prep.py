#!/usr/bin/env python3
"""
Assemble the no-motif family representatives (those matching NONE of their class's
8 STREME motifs, per FIMO) into a manifest (used by no_motif_umap_interactive.py).

Inputs (in this directory):
  reps_{cls}.clean.fasta     all family reps per class (what STREME/FIMO saw)
  fimo_{cls}/fimo.tsv         FIMO hits of the 8 class motifs vs those reps
Output (here):
  no_motif_manifest.csv       Entry, class, length (one row per no-motif rep)
"""
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
MA = HERE
CLASSES = ["DNA", "Protein", "RNA", "Metabolite"]


def read_fasta(path):
    seqs, sid, buf = {}, None, []
    for ln in open(path):
        ln = ln.rstrip()
        if ln.startswith(">"):
            if sid:
                seqs[sid] = "".join(buf)
            sid, buf = ln[1:].split()[0], []
        elif ln:
            buf.append(ln)
    if sid:
        seqs[sid] = "".join(buf)
    return seqs


def main():
    rows = []
    for cls in CLASSES:
        reps = read_fasta(MA / f"reps_{cls}.clean.fasta")
        fimo = pd.read_csv(MA / f"fimo_{cls}" / "fimo.tsv", sep="\t", comment="#")
        matched = set(fimo["sequence_name"])
        no_motif = [e for e in reps if e not in matched]
        for e in no_motif:
            rows.append({"Entry": e, "class": cls, "length": len(reps[e])})
        print(f"  {cls}: {len(no_motif)} no-motif reps")
    pd.DataFrame(rows).to_csv(HERE / "no_motif_manifest.csv", index=False)
    print(f"[done] {len(rows)} no-motif reps -> no_motif_manifest.csv")


if __name__ == "__main__":
    main()
