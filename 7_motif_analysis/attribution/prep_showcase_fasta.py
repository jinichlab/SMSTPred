#!/usr/bin/env python3
"""
Write a FASTA of selected SAM_swiss entries (sequences pulled from the labelled
CSV) for per-residue ProstT5 embedding + Integrated-Gradients attribution.

Default entry list = a curated single-label showcase, one (or two) per substrate
class, chosen so the IG map can be validated against a *known* catalytic motif:

  Q67472  DNA  (cytosine-5) MTase    -> Motif IV  "PC" catalytic dipeptide
  P04392  DNA  Dam, N6-adenine       -> "DPPY/NPPY" catalytic tetrad
  A2AA28  Protein  METTL23 (PRMT)    -> arginine-MTase double-E / THW
  Q9ZM61  RNA  RsmG (16S rRNA)       -> rRNA small-subunit MTase
  P9WEP8  Other  O-methyltransferase -> small-molecule O-MTase
"""
import argparse
from pathlib import Path

import pandas as pd

CSV = Path("/home/nroethler/Code/files/SAM_swiss_enzymes_substrate_types.csv")
SHOWCASE = ["Q67472", "P04392", "A2AA28", "Q9ZM61", "P9WEP8"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--entries", default=",".join(SHOWCASE),
                    help="comma-separated UniProt entries (default: curated showcase)")
    ap.add_argument("--csv", type=Path, default=CSV)
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "showcase.fasta")
    args = ap.parse_args()

    want = [e.strip() for e in args.entries.split(",") if e.strip()]
    df = pd.read_csv(args.csv, low_memory=False).drop_duplicates("Entry").set_index("Entry")

    n = 0
    with open(args.out, "w") as fh:
        for e in want:
            if e not in df.index:
                print(f"  [warn] {e} not in CSV, skipped")
                continue
            seq = str(df.loc[e, "Sequence"])
            name = str(df.loc[e, "Protein names"])[:60]
            fh.write(f">{e} {name}\n{seq}\n")
            n += 1
            print(f"  {e}  len={len(seq)}  {name}")
    print(f"[done] wrote {n} sequences -> {args.out}")


if __name__ == "__main__":
    main()
