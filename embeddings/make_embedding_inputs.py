#!/usr/bin/env python3
"""
What to embed: the FASTA files for the sequence embedders (ProtT5, ProstT5), written from the typed label tables of the newest pull.

  python make_embedding_inputs.py [--pull-dir DIR] [--release 2026_03] [--out <release>/embedding_inputs]

  sam_deployment.fasta   the SAM proteins with in_deployment (<= 2,000 aa) of labels_SAM.csv   -> embeddings_SAM_{ProtT5,ProstT5}.h5
  up_deployment.fasta    the UP enzymes with in_deployment of labels_UP.csv                    -> embeddings_UP_{ProtT5,ProstT5}.h5

The files hold every wanted protein, not only the new ones: the embedders compare each with what their output file already holds (incremental.py) and embed only what is new or
changed. The structure embedders (SaProt, ESM-IF) do not need a FASTA: they read the structure store of fetch_structures.py (the benchmark set, <= 1,000 aa).
Sequences are the pulled UniProt sequences (sam_mtases.tsv.gz / up_enzymes.tsv.gz of the same release), so a sequence that UniProt changed since the last pull is seen as changed.
"""
import argparse
import os
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
PULL = Path(os.environ.get("UNIPROT_PULL_DIR", REPO / "data" / "uniprot_pull"))
SETS = {"sam": ("labels_SAM.csv", "sam_mtases.tsv.gz"), "up": ("labels_UP.csv", "up_enzymes.tsv.gz")}


def newest_release(pull_dir, release=None):
    if release:
        return Path(pull_dir) / f"release_{release}"
    rel = sorted(p for p in Path(pull_dir).glob("release_*") if p.is_dir() and not p.name.startswith("_"))
    if not rel:
        sys.exit(f"no release_* folder in {pull_dir}")
    return rel[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pull-dir", default=str(PULL))
    ap.add_argument("--release")
    ap.add_argument("--out")
    ap.add_argument("--flag", default="in_deployment")
    a = ap.parse_args()
    rel = newest_release(a.pull_dir, a.release)
    out = Path(a.out) if a.out else rel / "embedding_inputs"
    out.mkdir(parents=True, exist_ok=True)
    for which, (lab_name, tsv) in SETS.items():
        lab = pd.read_csv(rel / lab_name, usecols=["Entry", a.flag], low_memory=False)
        want = lab.loc[lab[a.flag].astype(bool), "Entry"]
        seq = pd.read_csv(rel / tsv, sep="\t", usecols=["Entry", "Sequence"], low_memory=False).drop_duplicates("Entry").set_index("Entry")["Sequence"]
        seq = want.map(seq)
        bad = seq.isna() | (seq.astype(str).str.strip() == "")
        if bad.any():
            sys.exit(f"{int(bad.sum())} {which} proteins have no sequence, e.g. {list(want[bad][:5])}")
        path = out / f"{which}_{a.flag.replace('in_', '')}.fasta"
        with open(path, "w") as f:
            for e, s in zip(want, seq):
                f.write(f">{e}\n{s.strip().upper()}\n")
        lens = seq.str.len()
        print(f"{path}: {len(want)} sequences, length {lens.min()} / {int(lens.median())} / {lens.max()} (min / median / max)", flush=True)


if __name__ == "__main__":
    main()
