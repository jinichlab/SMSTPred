#!/usr/bin/env python3
"""
Stage 6 / motif analysis -- step 1: assemble the labeled SAM-MTase corpus.

Combines every SAM-MTase sequence with its substrate label:
  - training (7,405 SAM_swiss):   GROUND-TRUTH label (Metabolite/Protein/RNA/DNA cols)
  - deployment (tree of life):    PREDICTED label (call_* from predictions_output)

Each protein is assigned a single substrate (exactly one positive -> that class;
>1 -> "Mixed"; 0 -> "None"). Writes a master table (sequence + InterPro + class +
source + in_domain) and one FASTA per single-label class for downstream
de-redundification, domain enrichment, and STREME.

  python 1_build_corpus.py                      the archived inputs (data/labels_SAM_full.csv, data/reference_proteomes, 3_predictions/predictions_output); writes here
  python 1_build_corpus.py --release <folder>   a refresh release: <release>/labels_SAM.csv (in_deployment proteins, sequences + InterPro from sam_mtases.tsv.gz),
                                                <release>/predictions_all_proteomes.csv.gz (the ProtT5 ensemble's calls) + <release>/reference_proteomes/*.csv;
                                                writes <release>/motif_analysis/. (default release: $RELEASE_DIR, else the newest release_* of $UNIPROT_PULL_DIR)

Outputs (./, or <release>/motif_analysis/):
  corpus.csv                 Entry, Sequence, InterPro, class, source, in_domain
  corpus_all.fasta           all single-label sequences (header: >Entry|class)
  class_{RNA,DNA,Protein,Metabolite}.fasta
"""
import argparse
import glob
import os
import sys
from pathlib import Path

import pandas as pd

LABELS = ["Metabolite", "Protein", "RNA", "DNA"]
REPO = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("SMSTPRED_DATA", REPO / "data"))
TRAIN_CSV = str(DATA / "labels_SAM_full.csv")
REF_DIR = str(DATA / "reference_proteomes")      # per-organism UniProt proteome tables (not archived; ~30 GB)
PRED_DIR = str(REPO / "3_predictions" / "predictions_output")
PRED_SUFFIX = "_embeddings_predictions.csv"
OUT = Path(__file__).resolve().parent


def single_label(frame, cols):
    """Per row: the lone positive label, else 'Mixed' (>1) / 'None' (0)."""
    pos = frame[cols].apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    n = pos.sum(axis=1)
    lab = pos.idxmax(axis=1).map(dict(zip(cols, LABELS)))
    out = lab.where(n == 1, "Mixed")
    return out.where(n >= 1, "None")


def corpus_archived():
    """The archived inputs: SAM_swiss training table + the per-proteome prediction files of 3_predictions/predictions_output."""
    tr = pd.read_csv(TRAIN_CSV, low_memory=False)
    tr["class"] = single_label(tr, LABELS)
    train = tr[["Entry", "Sequence", "InterPro", "class"]].copy()
    train["source"] = "train"
    train["in_domain"] = 1
    print(f"training: {len(train)} | {train['class'].value_counts().to_dict()}", flush=True)

    csv_lookup = {Path(f).stem: f for f in glob.glob(f"{REF_DIR}/*.csv")}
    frames = []
    for pf in sorted(glob.glob(f"{PRED_DIR}/*{PRED_SUFFIX}")):
        stem = os.path.basename(pf)[:-len(PRED_SUFFIX)]
        ref = csv_lookup.get(stem)
        if ref is None:
            continue
        pred = pd.read_csv(pf, low_memory=False)
        callcols = [f"call_{l}" for l in LABELS]
        if not set(callcols).issubset(pred.columns):
            continue
        pred["class"] = single_label(pred, callcols)
        seqs = pd.read_csv(ref, usecols=["Entry", "Sequence", "InterPro"], low_memory=False)
        m = pred.merge(seqs, on="Entry", how="left")
        frames.append(m[["Entry", "Sequence", "InterPro", "class", "in_domain"]])
    return train, pd.concat(frames, ignore_index=True)


def corpus_release(rel):
    """A refresh release: training = the labelled SAM enzymes the deployment model covers (in_deployment, <= 2,000 aa), deployment = the proteome predictions."""
    lab = pd.read_csv(rel / "labels_SAM.csv", usecols=["Entry", *LABELS, "in_deployment"], low_memory=False)
    lab = lab[lab.in_deployment == True].copy()
    lab["class"] = single_label(lab, LABELS)
    seq = pd.read_csv(rel / "sam_mtases.tsv.gz", sep="\t", usecols=["Entry", "Sequence", "InterPro"], low_memory=False).drop_duplicates("Entry")
    train = lab[["Entry", "class"]].merge(seq, on="Entry", how="left")[["Entry", "Sequence", "InterPro", "class"]]
    train["source"] = "train"
    train["in_domain"] = 1
    print(f"training: {len(train)} | {train['class'].value_counts().to_dict()}", flush=True)

    pred = pd.read_csv(rel / "predictions_all_proteomes.csv.gz", low_memory=False)
    pred["class"] = single_label(pred, [f"call_{l}" for l in LABELS])
    ref = pd.concat([pd.read_csv(f, usecols=["Entry", "Sequence", "InterPro"], low_memory=False) for f in sorted((rel / "reference_proteomes").glob("*.csv"))],
                    ignore_index=True).drop_duplicates("Entry")
    dep = pred.merge(ref, on="Entry", how="left")[["Entry", "Sequence", "InterPro", "class", "in_domain"]]
    print(f"proteome predictions: {len(pred)}, with a sequence in reference_proteomes: {int(dep.Sequence.notna().sum())}", flush=True)
    return train, dep


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", nargs="?", const="", default=None, help="a refresh release folder (see above); without a value: $RELEASE_DIR or the newest release")
    a = ap.parse_args()
    out = OUT
    if a.release is not None:
        sys.path.insert(0, str(REPO / "1_clustering"))
        from mmseqs_clusters import newest_release
        rel = Path(a.release or os.environ.get("RELEASE_DIR") or newest_release(Path(os.environ.get("UNIPROT_PULL_DIR", REPO / "data" / "uniprot_pull")), os.environ.get("RELEASE")))
        out = rel / "motif_analysis"
        out.mkdir(exist_ok=True)
        train, dep = corpus_release(rel)
    else:
        train, dep = corpus_archived()
    dep["source"] = "deploy"
    print(f"deployment: {len(dep)} | {dep['class'].value_counts().to_dict()}", flush=True)

    # --- combine; prefer training rows on Entry collisions; require a sequence ---
    corpus = (pd.concat([train, dep], ignore_index=True)
              .dropna(subset=["Sequence"])
              .drop_duplicates("Entry", keep="first")
              .reset_index(drop=True))
    corpus["Sequence"] = corpus["Sequence"].str.replace(r"\s+", "", regex=True)
    corpus.to_csv(out / "corpus.csv", index=False)
    print(f"\ncorpus: {len(corpus)} unique proteins | "
          f"{corpus['class'].value_counts().to_dict()}", flush=True)

    # --- FASTAs (single-label classes only) ---
    sl = corpus[corpus["class"].isin(LABELS)]
    with open(out / "corpus_all.fasta", "w") as fa:
        for _, r in sl.iterrows():
            fa.write(f">{r.Entry}|{r['class']}\n{r.Sequence}\n")
    for lab in LABELS:
        sub = sl[sl["class"] == lab]
        with open(out / f"class_{lab}.fasta", "w") as fa:
            for _, r in sub.iterrows():
                fa.write(f">{r.Entry}\n{r.Sequence}\n")
        print(f"  class_{lab}.fasta: {len(sub)} sequences", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
