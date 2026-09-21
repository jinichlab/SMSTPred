#!/usr/bin/env python3
"""
Build the REAL per-protein motif-status dataset from FIMO run over ALL members
(not just the 50%-id representatives, and not homology-propagated).

For each substrate class C, FIMO scanned every class_C member against C's 8 STREME
motifs at p<1e-4 (fimo_all_C/fimo.tsv). A protein is "no-motif" if it matches NONE
of its own class's motifs.

Inputs (this dir):
  class_{DNA,RNA,Protein,Other}.fasta   all single-label members per class
  fimo_all_{cls}/fimo.tsv                FIMO hits of that class's 8 motifs vs members
Outputs:
  motif_status_all.csv        every protein: Entry, class, length, n_hits,
                              n_distinct_motifs, has_motif
  no_motif_real.csv           the no-motif subset only: Entry, class, length
"""
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
CLASSES = ["DNA", "RNA", "Protein", "Other"]


def read_fasta_lengths(path):
    lens, sid, n = {}, None, 0
    for ln in open(path):
        if ln.startswith(">"):
            if sid is not None:
                lens[sid] = n
            sid, n = ln[1:].split()[0].split("|")[0], 0
        else:
            n += len(ln.strip())
    if sid is not None:
        lens[sid] = n
    return lens


def main():
    rows = []
    for cls in CLASSES:
        lens = read_fasta_lengths(HERE / f"class_{cls}.fasta")
        fimo = pd.read_csv(HERE / f"fimo_all_{cls}" / "fimo.tsv", sep="\t", comment="#")
        fimo = fimo.dropna(subset=["sequence_name"])
        n_hits = fimo.groupby("sequence_name").size()
        n_distinct = fimo.groupby("sequence_name")["motif_id"].nunique()
        matched = set(fimo["sequence_name"])
        for e, L in lens.items():
            rows.append({"Entry": e, "class": cls, "length": L,
                         "n_hits": int(n_hits.get(e, 0)),
                         "n_distinct_motifs": int(n_distinct.get(e, 0)),
                         "has_motif": e in matched})
        nm = sum(1 for e in lens if e not in matched)
        print(f"  {cls}: {len(lens):>6} members | {len(lens)-nm:>6} with motif | "
              f"{nm:>6} no-motif ({100*nm/len(lens):.1f}%)")

    df = pd.DataFrame(rows)
    df.to_csv(HERE / "motif_status_all.csv", index=False)
    nomotif = df[~df["has_motif"]][["Entry", "class", "length"]]
    nomotif.to_csv(HERE / "no_motif_real.csv", index=False)
    tot, nm = len(df), int((~df["has_motif"]).sum())
    print(f"\nTOTAL: {tot:,} proteins | {tot-nm:,} carry their class motif | "
          f"{nm:,} no-motif ({100*nm/tot:.1f}%)")
    print(f"[wrote] motif_status_all.csv ({tot:,} rows) + no_motif_real.csv ({nm:,} rows)")


if __name__ == "__main__":
    main()
