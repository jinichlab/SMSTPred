#!/usr/bin/env python3
"""ProtT5-XL per-protein embeddings (1024-D) from a FASTA file.

  python embed_prott5.py --fasta proteins.fasta --out embeddings.h5

Output: HDF5 file, one dataset per FASTA id (the first word of the header; '/' and '.' become '_'), shape (1024,).

Incremental (incremental.py): an existing --out file is appended to, never rewritten. Only proteins that are not in it yet, or whose sequence changed since their vector was made,
are embedded; every vector carries the md5 of its sequence (attribute input_md5). For the archived files, which have no fingerprints, give --ledger (a FASTA / CSV of the sequences
they were made from): matching vectors are kept and stamped, the others embedded again. --overwrite starts the file again. The model is not even loaded when nothing is to be done.
  python embed_prott5.py --fasta wanted.fasta --out embeddings_SAM_ProtT5.h5 --ledger archived_sequences.csv

Two pipelines produced the archived ProtT5 files (see README.md):
  SAM and UP sets   --model Rostlab/prot_t5_xl_half_uniref50-enc --pool residue-mean --precision half     (float16 files; the default)
  tree-of-life set  --model Rostlab/prot_t5_xl_uniref50          --pool masked-mean  --precision full     (float32 files)
"""
import argparse
import re
import time
from pathlib import Path

import h5py
import numpy as np
import torch
from transformers import T5EncoderModel, T5Tokenizer

import incremental as inc


def id_norm(sid):
    return sid.replace("/", "_").replace(".", "_")


def read_fasta(path):
    seqs, sid = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            sid = id_norm(line[1:].split()[0])
            seqs[sid] = ""
        elif line and sid is not None:
            seqs[sid] += "".join(line.split()).upper().replace("-", "")
    return seqs


def batches(items, max_residues, max_seqs):
    """Longest first; a batch is closed when it would exceed max_residues (padded) or max_seqs."""
    items = sorted(items, key=lambda kv: -len(kv[1]))
    batch, longest = [], 0
    for sid, seq in items:
        if batch and (max(longest, len(seq)) * (len(batch) + 1) > max_residues or len(batch) >= max_seqs):
            yield batch
            batch, longest = [], 0
        batch.append((sid, seq))
        longest = max(longest, len(seq))
    if batch:
        yield batch


@torch.inference_mode()
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--model", default="Rostlab/prot_t5_xl_half_uniref50-enc")
    ap.add_argument("--pool", choices=["residue-mean", "masked-mean"], default="residue-mean",
                    help="residue-mean: mean over the L residue tokens; masked-mean: mean over every non-padding token, end-of-sequence token included")
    ap.add_argument("--precision", choices=["half", "full"], default="half",
                    help="half: fp16 model on a GPU and float16 output; full: fp32 model and float32 output (CPU always runs fp32)")
    ap.add_argument("--max-residues", type=int, default=4000, help="cap on batch size x longest sequence")
    ap.add_argument("--max-seqs", type=int, default=100)
    ap.add_argument("--ledger", nargs="+", help="what the unstamped vectors already in --out were made from (FASTA, or CSV / TSV with Entry and Sequence)")
    ap.add_argument("--trust-existing", action="store_true", help="keep unstamped vectors by accession alone, without a ledger")
    ap.add_argument("--overwrite", action="store_true", help="start --out again instead of appending")
    args = ap.parse_args()

    seqs = read_fasta(args.fasta)
    wanted = {sid: inc.md5_text(seq) for sid, seq in seqs.items()}
    if args.overwrite or not args.out.exists():
        pl = {"todo": {sid: "new" for sid in wanted}, "keep": [], "stamp": {}, "extra": []}
    else:
        pl = inc.plan(wanted, inc.stored_fingerprints(args.out), inc.read_ledger(args.ledger, norm=id_norm), args.trust_existing)
    inc.report("ProtT5", pl, len(seqs))
    if not pl["todo"]:
        inc.apply_stamps(args.out, pl["stamp"])
        print("[done] nothing to embed", flush=True)
        return

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    half = args.precision == "half" and device.type == "cuda"
    print(f"[load] {args.model} on {device} ({'fp16' if half else 'fp32'}), pooling {args.pool}", flush=True)
    tok = T5Tokenizer.from_pretrained(args.model, do_lower_case=False)
    model = T5EncoderModel.from_pretrained(args.model, torch_dtype=torch.float16 if half else torch.float32).to(device).eval()

    todo = [(sid, seqs[sid]) for sid in pl["todo"]]
    print(f"[data] {len(todo)} sequences to embed", flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    t0, done = time.time(), 0
    with h5py.File(args.out, "w" if args.overwrite else "a") as h5:
        inc.check_config(h5, {"model": args.model, "pool": args.pool, "precision": "half" if half else "full"})
        for batch in batches(todo, args.max_residues, args.max_seqs):
            texts = [" ".join(re.sub(r"[UZOB]", "X", s)) for _, s in batch]
            enc = tok(texts, add_special_tokens=True, padding="longest", return_tensors="pt").to(device)
            h = model(input_ids=enc.input_ids, attention_mask=enc.attention_mask).last_hidden_state
            for i, (sid, seq) in enumerate(batch):
                if args.pool == "residue-mean":
                    v = h[i, :len(seq)].mean(dim=0)
                else:
                    m = enc.attention_mask[i].bool()
                    v = h[i][m].mean(dim=0)
                inc.replace_slot(h5, sid)
                ds = h5.create_dataset(sid, data=v.float().cpu().numpy().astype(np.float16 if half else np.float32))
                ds.attrs[inc.ATTR] = wanted[sid]
            done += len(batch)
            print(f"[progress] {done}/{len(todo)} ({done / max(time.time() - t0, 1e-9):.1f} seq/s)", flush=True)
    inc.apply_stamps(args.out, pl["stamp"])
    print(f"[done] {args.out}: embedded {len(todo)}, kept {len(pl['keep'])}", flush=True)


if __name__ == "__main__":
    main()
