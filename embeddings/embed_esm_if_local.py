#!/usr/bin/env python3
"""
ESM-IF embeddings from LOCAL structure files, appended to an existing H5.

Companion to embed_esm_if.py (which streams structures from AlphaFold DB). Use
this for proteins absent from AFDB whose structures were predicted locally
(e.g. AF2 / localcolabfold). For each protein, the best-model PDB in
--struct-dir is loaded, ESM-IF encodes chain A and mean-pools to 512-D, and the
embedding is written into --out keyed by UniProt Entry. The Entry is taken as
the filename prefix before the first '_' (ColabFold names files
"<Entry>_..._rank_001_..._.pdb"; a plain "<Entry>.pdb" works too).

Incremental (incremental.py): --out is appended to. Only structures that are not in it yet, or whose PDB file changed since their vector was made (a re-predicted structure, a new
AlphaFold DB model), are embedded; every vector carries the md5 of its PDB file (attribute input_md5). For the archived file, which has no fingerprints, give --ledger <folder of the
archived <Entry>.pdb files>: identical files are kept and stamped, the others embedded again. This is the ESM-IF embedder of the refresh: it reads the structure store of
fetch_structures.py, so ESM-IF and SaProt see the same structures. (embed_esm_if.py, which streams AlphaFold DB, made the archived file; it keeps no structure and checks no sequence.)

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
import incremental as inc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--struct-dir', required=True,
                    help='directory of predicted PDBs (e.g. ColabFold output)')
    ap.add_argument('--out', required=True, help='ESM-IF H5 to append into')
    ap.add_argument('--pattern', default='*rank_001*.pdb',
                    help='glob for the best-model file per protein')
    ap.add_argument('--chain', default='A')
    ap.add_argument('--device', choices=['cuda', 'cpu'], default='cuda')
    ap.add_argument('--ledger', nargs='+', help='folder(s) of the <Entry>.pdb files the unstamped vectors already in --out were made from')
    ap.add_argument('--trust-existing', action='store_true', help='keep unstamped vectors by accession alone, without a ledger')
    ap.add_argument('--overwrite', action='store_true', help='start --out again instead of appending')
    args = ap.parse_args()

    # Entry -> best-model PDB (Entry = filename prefix before first '_').
    ent2pdb = {}
    for f in sorted(glob.glob(os.path.join(args.struct_dir, args.pattern))):
        entry = os.path.splitext(os.path.basename(f))[0].split('_')[0]      # <Entry>.pdb or ColabFold's <Entry>_..._rank_001_....pdb
        ent2pdb.setdefault(entry, f)
    print(f"found {len(ent2pdb)} proteins matching '{args.pattern}' in {args.struct_dir}",
          flush=True)

    wanted = {e: inc.md5_file(p) for e, p in ent2pdb.items()}
    if args.overwrite or not Path(args.out).exists():
        pl = {"todo": {e: "new" for e in wanted}, "keep": [], "stamp": {}, "extra": []}
    else:
        pl = inc.plan(wanted, inc.stored_fingerprints(args.out), inc.read_ledger(args.ledger), args.trust_existing)
    inc.report("ESM-IF", pl, len(wanted))
    todo = {e: ent2pdb[e] for e in pl["todo"]}
    if not todo:
        inc.apply_stamps(args.out, pl["stamp"])
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
    with h5py.File(args.out, 'w' if args.overwrite else 'a') as h5:
        inc.check_config(h5, {'model': 'esm_if1_gvp4_t16_142M_UR50', 'chain': args.chain})
        for entry, pdb in todo.items():
            try:
                coords, _ = esm.inverse_folding.util.load_coords(pdb, args.chain)
                device = next(model.parameters()).device
                rep = encode_structure(model, alphabet, coords, device)
                emb = rep.mean(dim=0).detach().cpu().numpy().astype(np.float32)
                inc.replace_slot(h5, entry)
                h5.create_dataset(entry, data=emb, compression='gzip').attrs[inc.ATTR] = wanted[entry]
                n_ok += 1
                print(f"  {entry}: L={len(coords)} dim={emb.shape[0]}", flush=True)
            except Exception as e:
                n_fail += 1
                print(f"  ! {entry}: {e}", flush=True)

    inc.apply_stamps(args.out, pl["stamp"])
    print(f"\nembedded {n_ok}, failed {n_fail}, kept {len(pl['keep'])}", flush=True)
    with h5py.File(args.out, 'r') as f:
        print(f"{args.out}: {len(f.keys())} embeddings total", flush=True)


if __name__ == '__main__':
    main()
