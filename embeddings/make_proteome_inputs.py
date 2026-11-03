#!/usr/bin/env python3
"""
What to embed for the reference proteomes: the methyltransferase-domain proteins (InterPro IPR029063) up to the model's length limit.

  python make_proteome_inputs.py [--pull-dir DIR] [--release 2026_03] [--max-len 2000] [--out <release>/embedding_inputs]

Reads <release>/reference_proteomes/<proteome>.csv (every protein of a proteome up to 2,000 aa) and <release>/reference_proteomes_long/<proteome>.csv (longer ones). Writes
  proteomes_ipr029063.fasta         the proteins with IPR029063 and Length <= --max-len (default 2,000 aa: the longest training protein of the deployment ensemble is 1,848 aa), one record per Entry
  proteome_members.csv              proteome (file stem), Entry, Length of every protein in the fasta -> split_proteome_embeddings.py, 3_predictions
  proteomes_not_predicted.csv       proteins with IPR029063 LONGER than --max-len: outside what the model has seen; they get no prediction, and are listed here
"""
import argparse
import glob
import os
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
PULL = Path(os.environ.get("UNIPROT_PULL_DIR", HERE.parent / "data" / "uniprot_pull"))
DOMAIN = "IPR029063"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pull-dir", default=str(PULL))
    ap.add_argument("--release")
    ap.add_argument("--max-len", type=int, default=2000)
    ap.add_argument("--out")
    a = ap.parse_args()
    if os.environ.get("RELEASE_DIR"):
        rel = Path(os.environ["RELEASE_DIR"])
    elif a.release:
        rel = Path(a.pull_dir) / f"release_{a.release}"
    else:
        cand = sorted(p for p in Path(a.pull_dir).glob("release_*") if p.is_dir() and not p.name.startswith("_"))
        if not cand:
            sys.exit(f"no release_* folder in {a.pull_dir}")
        rel = cand[-1]
    out = Path(a.out) if a.out else rel / "embedding_inputs"
    out.mkdir(parents=True, exist_ok=True)
    keep, skipped = [], []
    for d in ("reference_proteomes", "reference_proteomes_long"):
        files = sorted(glob.glob(str(rel / d / "*.csv")))
        if not files:
            print(f"WARNING: no csv files in {rel / d}", flush=True)
        for f in files:
            t = pd.read_csv(f, usecols=["Entry", "Length", "Sequence", "InterPro"], low_memory=False)
            t = t[t.InterPro.fillna("").str.contains(DOMAIN)].assign(proteome=Path(f).stem)
            keep.append(t[t.Length <= a.max_len])
            skipped.append(t[t.Length > a.max_len][["proteome", "Entry", "Length"]])
    k = pd.concat(keep, ignore_index=True)
    k = k[k.Sequence.notna() & (k.Sequence.str.strip() != "")]
    members = k[["proteome", "Entry", "Length"]]
    uniq = k.drop_duplicates("Entry")
    with open(out / "proteomes_ipr029063.fasta", "w") as fh:
        for e, s in zip(uniq.Entry, uniq.Sequence):
            fh.write(f">{e}\n{s.strip().upper()}\n")
    members.to_csv(out / "proteome_members.csv", index=False)
    pd.concat(skipped).to_csv(out / "proteomes_not_predicted.csv", index=False)
    print(f"{len(members)} proteome entries ({len(uniq)} distinct proteins, {members.proteome.nunique()} proteomes) with {DOMAIN} and length <= {a.max_len} -> proteomes_ipr029063.fasta; "
          f"length {int(members.Length.min())} / {int(members.Length.median())} / {int(members.Length.max())} (min / median / max)")
    print(f"{sum(len(s) for s in skipped)} longer than {a.max_len} aa are not predicted -> proteomes_not_predicted.csv", flush=True)


if __name__ == "__main__":
    main()
