#!/usr/bin/env python3
"""
ESM-IF embeddings from LOCAL structure files, appended to an existing H5.

Companion to embed_esm_if.py (which streams structures from AlphaFold DB). Use
this for proteins absent from AFDB whose structures were predicted locally
(e.g. AF2 / localcolabfold). For each protein, the best-model PDB in
--struct-dir is loaded, ESM-IF encodes chain A and mean-pools to 512-D, and the
embedding is written into --out keyed by UniProt Entry. The Entry is taken as
the filename prefix before the first '_' (ColabFold names files
"<Entry>_..._rank_001_..._.pdb"; a plain "<Entry>.pdb" works too). Resumable: entries already in --out are
skipped, so it safely appends to an AFDB-derived H5.

Example (the 26 ColabFold structures of the SAM set, or any folder of <Entry>.pdb files):
  python embed_esm_if_local.py --struct-dir structures/ --pattern '*.pdb' \
      --out embeddings_SAM_ESM-IF.h5 --device cuda
"""

import argparse
import glob
import os
from pathlib import Path

import h5py
import numpy as np
import torch

import esm
import esm.inverse_folding

from embed_esm_if import encode_structure  # reuse the device-safe encoder


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--struct-dir', required=True,
                    help='directory of predicted PDBs (e.g. ColabFold output)')
    ap.add_argument('--out', required=True, help='ESM-IF H5 to append into')
    ap.add_argument('--pattern', default='*rank_001*.pdb',
                    help='glob for the best-model file per protein')
    ap.add_argument('--chain', default='A')
    ap.add_argument('--device', choices=['cuda', 'cpu'], default='cuda')
    args = ap.parse_args()

    # Entry -> best-model PDB (Entry = filename prefix before first '_').
    ent2pdb = {}
    for f in sorted(glob.glob(os.path.join(args.struct_dir, args.pattern))):
        entry = os.path.splitext(os.path.basename(f))[0].split('_')[0]      # <Entry>.pdb or ColabFold's <Entry>_..._rank_001_....pdb
        ent2pdb.setdefault(entry, f)
    print(f"found {len(ent2pdb)} proteins matching '{args.pattern}' in {args.struct_dir}",
          flush=True)

    done = set()
    if Path(args.out).exists():
        with h5py.File(args.out, 'r') as f:
            done = set(f.keys())
    todo = {e: p for e, p in ent2pdb.items() if e not in done}
    print(f"{len(done)} already in {args.out}; {len(todo)} to embed", flush=True)
    if not todo:
        print("nothing to do", flush=True)
        return

    print("loading ESM-IF1...", flush=True)
    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    model = model.eval()
    if args.device == 'cuda' and torch.cuda.is_available():
        model = model.cuda()
        print(f"  cuda: {torch.cuda.get_device_name(0)}", flush=True)
    else:
        print(f"  device: cpu (cuda_available={torch.cuda.is_available()})", flush=True)

    n_ok = n_fail = 0
    with h5py.File(args.out, 'a') as h5:
        for entry, pdb in todo.items():
            try:
                coords, _ = esm.inverse_folding.util.load_coords(pdb, args.chain)
                device = next(model.parameters()).device
                rep = encode_structure(model, alphabet, coords, device)
                emb = rep.mean(dim=0).detach().cpu().numpy().astype(np.float32)
                h5.create_dataset(entry, data=emb, compression='gzip')
                n_ok += 1
                print(f"  {entry}: L={len(coords)} dim={emb.shape[0]}", flush=True)
            except Exception as e:
                n_fail += 1
                print(f"  ! {entry}: {e}", flush=True)

    print(f"\nembedded {n_ok}, failed {n_fail}", flush=True)
    with h5py.File(args.out, 'r') as f:
        print(f"{args.out}: {len(f.keys())} embeddings total", flush=True)


if __name__ == '__main__':
    main()
