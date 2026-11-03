#!/usr/bin/env python3
"""One table of every proteome prediction of a refresh release.

  python combine_predictions.py [--release <folder>]      (default: $RELEASE_DIR, else the newest release_* of $UNIPROT_PULL_DIR)

Reads <release>/predictions/<proteome>_embeddings_predictions.csv (Stage 3, one file per proteome) and <release>/embedding_inputs/proteome_members.csv
(proteome, Entry, Length) and writes <release>/predictions_all_proteomes.csv.gz: the columns of the prediction files plus `proteome` (the file stem) and `Length`.
This is the input of 4_phylo_trends/reference_tree.py --release and of the Stage 7 corpus (7_motif_analysis/1_build_corpus.py --release).
"""
import argparse
import os
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
SUFFIX = "_embeddings_predictions.csv"


def combine(rel):
    files = sorted((rel / "predictions").glob(f"*{SUFFIX}"))
    if not files:
        sys.exit(f"no prediction files in {rel / 'predictions'}: run the proteome embedding and Stage 3 first")
    pred = pd.concat([pd.read_csv(f).assign(proteome=f.name[:-len(SUFFIX)]) for f in files], ignore_index=True)
    members = pd.read_csv(rel / "embedding_inputs" / "proteome_members.csv", usecols=["proteome", "Entry", "Length"])
    out = pred.merge(members, on=["proteome", "Entry"], how="left")
    assert len(out) == len(pred), "an Entry occurs twice in one proteome"
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release")
    a = ap.parse_args()
    sys.path.insert(0, str(REPO / "1_clustering"))
    from mmseqs_clusters import newest_release
    rel = Path(a.release or os.environ.get("RELEASE_DIR") or newest_release(Path(os.environ.get("UNIPROT_PULL_DIR", REPO / "data" / "uniprot_pull")), os.environ.get("RELEASE")))
    out = combine(rel)
    out.to_csv(rel / "predictions_all_proteomes.csv.gz", index=False)
    print(f"{len(out)} proteins from {out.proteome.nunique()} proteomes -> {rel}/predictions_all_proteomes.csv.gz ({int(out.Length.isna().sum())} without a length)")


if __name__ == "__main__":
    main()
