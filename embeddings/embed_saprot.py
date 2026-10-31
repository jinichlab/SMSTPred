#!/usr/bin/env python3
"""SaProt (650M, AF2) structure-aware embeddings for SAM_swiss.

For each protein: foldseek tokenises the predicted structure into 3Di, the
amino acid + 3Di are fused per residue into SaProt's structure-aware (SA)
sequence (uppercase AA + lowercase 3Di), and SaProt's encoder mean-pools the
residue tokens to a 1280-D vector. Unlike our ProstT5 embeddings (sequence in,
via <AA2fold>), SaProt is *conditioned on the actual predicted fold*.

Output mirrors the other embedding H5s: one float32 dataset per Entry. Needs the
structures from fetch_structures.py and the foldseek binary. Resumable.

Env: environment_transformers.yml (transformers + CUDA torch). Example:
  python embed_saprot.py --struct-dir data/sam_structures --out embeddings_SAM_SaProt.h5
"""
import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import h5py
import numpy as np
import torch
from transformers import EsmTokenizer, EsmForMaskedLM

MAX_TOKENS = 1024  # ESM-2 position limit; SA token == 1 residue


def get_sa_seq(foldseek, pdb):
    """Run foldseek 3Di descriptor; return the fused SA sequence (first chain).

    structureto3didescriptor emits one row per chain: name<TAB>AA<TAB>3Di<TAB>...
    AF2/AFDB models are single-chain, so we take the first data row."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "out")
        subprocess.run([foldseek, "structureto3didescriptor", "-v", "0",
                        "--threads", "1", "--chain-name-mode", "1", pdb, out],
                       check=True, capture_output=True)
        with open(out) as f:
            line = f.readline().rstrip("\n")
    if not line:
        return ""
    parts = line.split("\t")
    aa, tdi = parts[1], parts[2]
    return "".join(a + b.lower() for a, b in zip(aa, tdi))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--struct-dir", required=True, help="folder of <Entry>.pdb files (data/sam_structures after `make data`)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--model", default="westlake-repl/SaProt_650M_AF2")
    ap.add_argument("--foldseek",
                    default=os.environ.get("FOLDSEEK_BIN") or shutil.which("foldseek")
                    or str(Path(__file__).resolve().parents[1] / "tools" / "foldseek" / "bin" / "foldseek"),
                    help="Foldseek binary ($FOLDSEEK_BIN, else foldseek on PATH, else tools/foldseek/bin/foldseek from `make tools`)")
    ap.add_argument("--pattern", default="*.pdb")
    args = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"device={device}; loading {args.model}", flush=True)
    tok = EsmTokenizer.from_pretrained(args.model)
    model = EsmForMaskedLM.from_pretrained(args.model).to(device).eval()

    pdbs = sorted(Path(args.struct_dir).glob(args.pattern))
    done = set()
    if Path(args.out).exists():
        with h5py.File(args.out, "r") as f:
            done = set(f.keys())
    print(f"{len(pdbs)} structures; {len(done)} already embedded", flush=True)

    n_ok = n_fail = 0
    with h5py.File(args.out, "a") as h5:
        for i, pdb in enumerate(pdbs, 1):
            entry = pdb.stem
            if entry in done:
                continue
            try:
                sa = get_sa_seq(args.foldseek, str(pdb))
                if not sa:
                    raise ValueError("empty SA sequence")
                ids = tok(sa, return_tensors="pt", truncation=True,
                          max_length=MAX_TOKENS).to(device)
                with torch.no_grad():
                    h = model.esm(**ids).last_hidden_state[0, 1:-1]  # strip CLS/EOS
                emb = h.mean(0).cpu().numpy().astype(np.float32)
                h5.create_dataset(entry, data=emb, compression="gzip")
                n_ok += 1
            except Exception as e:
                n_fail += 1
                print(f"  ! {entry}: {e}", flush=True)
            if i % 500 == 0:
                print(f"  {i}/{len(pdbs)} (ok={n_ok} fail={n_fail})", flush=True)

    print(f"done: {n_ok} embedded, {n_fail} failed", flush=True)
    with h5py.File(args.out, "r") as f:
        ks = list(f.keys())
        dim = f[ks[0]].shape[0] if ks else 0
        print(f"{args.out}: {len(ks)} embeddings, dim={dim}", flush=True)


if __name__ == "__main__":
    main()
