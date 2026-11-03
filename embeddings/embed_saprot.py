#!/usr/bin/env python3
"""SaProt (650M, AF2) structure-aware embeddings for SAM_swiss.

For each protein: foldseek tokenises the predicted structure into 3Di, the
amino acid + 3Di are fused per residue into SaProt's structure-aware (SA)
sequence (uppercase AA + lowercase 3Di), and SaProt's encoder mean-pools the
residue tokens to a 1280-D vector. Unlike our ProstT5 embeddings (sequence in,
via <AA2fold>), SaProt is *conditioned on the actual predicted fold*.

Output mirrors the other embedding H5s: one float32 dataset per Entry. Needs the
structures from fetch_structures.py and the foldseek binary.

Incremental (incremental.py): an existing --out file is appended to. Only structures that are not in it yet, or whose PDB file changed since their vector was made (a re-predicted
structure, a new AlphaFold DB model), are embedded; every vector carries the md5 of its PDB file (attribute input_md5). For the archived file, which has no fingerprints, give --ledger
<folder of the archived <Entry>.pdb files>: identical files are kept and stamped, the others embedded again. The model is not loaded when nothing is to be done.

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

import incremental as inc

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
    ap.add_argument("--ledger", nargs="+", help="folder(s) of the <Entry>.pdb files the unstamped vectors already in --out were made from")
    ap.add_argument("--trust-existing", action="store_true", help="keep unstamped vectors by accession alone, without a ledger")
    ap.add_argument("--overwrite", action="store_true", help="start --out again instead of appending")
    args = ap.parse_args()

    pdbs = sorted(Path(args.struct_dir).glob(args.pattern))
    by_entry = {pdb.stem: pdb for pdb in pdbs}
    wanted = {e: inc.md5_file(f) for e, f in by_entry.items()}
    if args.overwrite or not Path(args.out).exists():
        pl = {"todo": {e: "new" for e in wanted}, "keep": [], "stamp": {}, "extra": []}
    else:
        pl = inc.plan(wanted, inc.stored_fingerprints(args.out), inc.read_ledger(args.ledger), args.trust_existing)
    inc.report("SaProt", pl, len(wanted))
    if not pl["todo"]:
        inc.apply_stamps(args.out, pl["stamp"])
        print("done: nothing to embed", flush=True)
        return

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"device={device}; loading {args.model}", flush=True)
    tok = EsmTokenizer.from_pretrained(args.model)
    model = EsmForMaskedLM.from_pretrained(args.model).to(device).eval()

    n_ok = n_fail = 0
    with h5py.File(args.out, "w" if args.overwrite else "a") as h5:
        inc.check_config(h5, {"model": args.model})
        for i, entry in enumerate(pl["todo"], 1):
            pdb = by_entry[entry]
            try:
                sa = get_sa_seq(args.foldseek, str(pdb))
                if not sa:
                    raise ValueError("empty SA sequence")
                ids = tok(sa, return_tensors="pt", truncation=True,
                          max_length=MAX_TOKENS).to(device)
                with torch.no_grad():
                    h = model.esm(**ids).last_hidden_state[0, 1:-1]  # strip CLS/EOS
                emb = h.mean(0).cpu().numpy().astype(np.float32)
                inc.replace_slot(h5, entry)
                h5.create_dataset(entry, data=emb, compression="gzip").attrs[inc.ATTR] = wanted[entry]
                n_ok += 1
            except Exception as e:
                n_fail += 1
                print(f"  ! {entry}: {e}", flush=True)
            if i % 500 == 0:
                print(f"  {i}/{len(pl['todo'])} (ok={n_ok} fail={n_fail})", flush=True)

    inc.apply_stamps(args.out, pl["stamp"])
    print(f"done: {n_ok} embedded, {n_fail} failed, {len(pl['keep'])} kept", flush=True)
    with h5py.File(args.out, "r") as f:
        ks = list(f.keys())
        dim = f[ks[0]].shape[0] if ks else 0
        print(f"{args.out}: {len(ks)} embeddings, dim={dim}", flush=True)


if __name__ == "__main__":
    main()
