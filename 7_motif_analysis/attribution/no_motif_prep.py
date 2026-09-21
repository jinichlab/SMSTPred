#!/usr/bin/env python3
"""
Assemble the no-motif family representatives (those matching NONE of their class's
8 STREME motifs, per FIMO) into a FASTA + manifest for per-residue IG attribution.

Inputs (in ../ = motif_analysis/):
  reps_{cls}.clean.fasta     all family reps per class (what STREME/FIMO saw)
  fimo_{cls}/fimo.tsv         FIMO hits of the 8 class motifs vs those reps
Outputs (here):
  no_motif.fasta             one record per no-motif rep, header ">Entry class"
  no_motif_manifest.csv       Entry, class, length
"""
from pathlib import Path
import pandas as pd

MA = Path(__file__).resolve().parent.parent      # motif_analysis/
HERE = Path(__file__).resolve().parent
CLASSES = ["DNA", "Protein", "RNA", "Other"]


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
    with open(HERE / "no_motif.fasta", "w") as out:
        for cls in CLASSES:
            reps = read_fasta(MA / f"reps_{cls}.clean.fasta")
            fimo = pd.read_csv(MA / f"fimo_{cls}" / "fimo.tsv", sep="\t", comment="#")
            matched = set(fimo["sequence_name"])
            no_motif = [e for e in reps if e not in matched]
            for e in no_motif:
                out.write(f">{e} {cls}\n{reps[e]}\n")
                rows.append({"Entry": e, "class": cls, "length": len(reps[e])})
            print(f"  {cls}: {len(no_motif)} no-motif reps")
    pd.DataFrame(rows).to_csv(HERE / "no_motif_manifest.csv", index=False)
    print(f"[done] {len(rows)} no-motif reps -> no_motif.fasta + no_motif_manifest.csv")


if __name__ == "__main__":
    main()
