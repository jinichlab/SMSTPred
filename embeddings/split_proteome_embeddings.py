#!/usr/bin/env python3
"""
Split the combined proteome embedding file into one H5 per proteome, as the archived tree-of-life files are (<proteome>_embeddings.h5, Entry -> vector), which 3_predictions/predict.py scores file by file.

  python split_proteome_embeddings.py [--pull-dir DIR] [--release 2026_03]
Reads <release>/embeddings/embeddings_proteomes_ProtT5.h5 and <release>/embedding_inputs/proteome_members.csv; writes <release>/embeddings/proteomes_ProtT5/<proteome>_embeddings.h5.
An Entry listed in several proteomes is written to each. Stops if a wanted protein has no vector.
"""
import argparse
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PULL = Path(os.environ.get("UNIPROT_PULL_DIR", HERE.parent / "data" / "uniprot_pull"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pull-dir", default=str(PULL))
    ap.add_argument("--release")
    a = ap.parse_args()
    if os.environ.get("RELEASE_DIR"):
        rel = Path(os.environ["RELEASE_DIR"])
    elif a.release:
        rel = Path(a.pull_dir) / f"release_{a.release}"
    else:
        rel = sorted(p for p in Path(a.pull_dir).glob("release_*") if p.is_dir() and not p.name.startswith("_"))[-1]
    members = pd.read_csv(rel / "embedding_inputs" / "proteome_members.csv")
    out = rel / "embeddings" / "proteomes_ProtT5"
    out.mkdir(parents=True, exist_ok=True)
    n = 0
    with h5py.File(rel / "embeddings" / "embeddings_proteomes_ProtT5.h5", "r") as src:
        have = set(src.keys())
        miss = sorted(set(members.Entry.str.replace("/", "_").str.replace(".", "_")) - have)
        if miss:
            sys.exit(f"{len(miss)} proteome proteins have no vector, e.g. {miss[:5]}: rerun embeddings/run_embeddings.sh (EMBED_WHAT=proteomes)")
        for stem, g in members.groupby("proteome"):
            with h5py.File(out / f"{stem}_embeddings.h5", "w") as dst:
                for e in g.Entry:
                    dst.create_dataset(e, data=np.asarray(src[e.replace("/", "_").replace(".", "_")][:], dtype=np.float32))
                    n += 1
    print(f"{n} vectors written into {members.proteome.nunique()} files in {out}", flush=True)


if __name__ == "__main__":
    main()
