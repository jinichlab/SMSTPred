#!/usr/bin/env python3
"""Re-embed a small random sample of the SAM training proteins with the scripts of this folder and compare with the archived vectors.

The sample is 20 proteins (seeded) that have a Swiss-Prot sequence of at most 1,000 residues, a structure in data/sam_structures and a vector in
all four archived SAM embedding files. Each model is run exactly as a user would run it (the embed_*.py scripts, as subprocesses) and the new vector
is compared with the one in the Zenodo archive: cosine similarity and maximum absolute difference per protein.

  python spotcheck_embeddings.py --models prott5 prostt5 saprot     # environment_transformers.yml, needs Foldseek for SaProt
  python spotcheck_embeddings.py --models esmif                     # environment_esmif.yml

Writes spotcheck_output/spotcheck_report.csv (one row per model and protein) and spotcheck_output/summary_<models>.txt. Exit status 1 if any model's minimum
cosine similarity is below --min-cosine.
"""
import argparse
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = Path(os.environ.get("SMSTPRED_DATA", HERE.parent / "data"))
OUT = HERE / "spotcheck_output"

ARCHIVE = {"prott5": "embeddings_SAM_ProtT5.h5", "prostt5": "embeddings_SAM_ProstT5.h5",
           "saprot": "embeddings_SAM_SaProt.h5", "esmif": "embeddings_SAM_ESM-IF.h5"}


def archived_vector(model, acc):
    with h5py.File(DATA / ARCHIVE[model]) as h:
        x = h[acc]
        return np.asarray(x["mean"] if isinstance(x, h5py.Group) else x, dtype=np.float64)


def sample(n, seed):
    lab = pd.read_csv(DATA / "labels_SAM_full.csv", usecols=["Entry", "Length", "Sequence"])
    lab = lab[lab.Length <= 1000]
    have = None
    for f in ARCHIVE.values():
        with h5py.File(DATA / f) as h:
            k = set(h.keys())
        have = k if have is None else have & k
    ok = [e for e in sorted(lab.Entry) if e in have and (DATA / "sam_structures" / f"{e}.pdb").exists()]
    pick = sorted(random.Random(seed).sample(ok, n))
    return lab.set_index("Entry").loc[pick].Sequence.to_dict()


def run(cmd):
    print("$", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True, cwd=HERE)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", nargs="+", choices=list(ARCHIVE), default=list(ARCHIVE))
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min-cosine", type=float, default=0.99)
    ap.add_argument("--foldseek", default=None, help="Foldseek binary for SaProt (default: see embed_saprot.py)")
    args = ap.parse_args()

    seqs = sample(args.n, args.seed)
    work = OUT / "work"
    shutil.rmtree(work, ignore_errors=True)
    (work / "structures").mkdir(parents=True)
    with open(work / "sample.fasta", "w") as f:
        for acc, s in seqs.items():
            f.write(f">{acc}\n{s}\n")
    for acc in seqs:
        os.symlink((DATA / "sam_structures" / f"{acc}.pdb").resolve(), work / "structures" / f"{acc}.pdb")
    print(f"{len(seqs)} proteins: {', '.join(seqs)}", flush=True)

    py = sys.executable
    new = {}
    for m in args.models:
        out = work / f"{m}.h5"
        if m == "prott5":
            run([py, "embed_prott5.py", "--fasta", work / "sample.fasta", "--out", out])
        elif m == "prostt5":
            run([py, "embed_prostt5.py", "--fasta", work / "sample.fasta", "--out", out, "--no-per-residue"])
        elif m == "saprot":
            run([py, "embed_saprot.py", "--struct-dir", work / "structures", "--out", out] + (["--foldseek", args.foldseek] if args.foldseek else []))
        elif m == "esmif":
            run([py, "embed_esm_if_local.py", "--struct-dir", work / "structures", "--pattern", "*.pdb", "--out", out])
        new[m] = out

    rows = []
    for m, path in new.items():
        with h5py.File(path) as h:
            for acc in seqs:
                x = h[acc]
                v = np.asarray(x["mean"] if isinstance(x, h5py.Group) else x, dtype=np.float64)
                a = archived_vector(m, acc)
                rows.append(dict(model=m, entry=acc, length=len(seqs[acc]), cosine=float(v @ a / (np.linalg.norm(v) * np.linalg.norm(a))),
                                 max_abs_diff=float(np.abs(v - a).max()), rel_norm_diff=float(np.linalg.norm(v - a) / np.linalg.norm(a))))
    df = pd.DataFrame(rows)
    tag = "_".join(args.models)
    df.to_csv(OUT / f"spotcheck_report_{tag}.csv", index=False)
    lines = [f"{len(seqs)} proteins, seed {args.seed}; cosine similarity of the re-computed vector with the archived one"]
    ok = True
    for m, g in df.groupby("model", sort=False):
        good = g.cosine.min() >= args.min_cosine
        ok &= bool(good)
        lines.append(f"{m:8s} min cosine {g.cosine.min():.6f}  median {g.cosine.median():.6f}  max |diff| {g.max_abs_diff.max():.2e}  "
                     f"relative norm difference (median) {g.rel_norm_diff.median():.2e}   {'PASS' if good else 'FAIL'} (threshold {args.min_cosine})")
    (OUT / f"summary_{tag}.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    shutil.rmtree(work, ignore_errors=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
