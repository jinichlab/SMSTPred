"""
Compute ProstT5 per-residue and per-protein embeddings.

Run in an environment with transformers, sentencepiece and h5py (environment_transformers.yml):
    python embed_prostt5.py \
        --fasta input.fasta --out embeddings.h5

Input: FASTA of amino-acid sequences.
Output: HDF5 file with groups
    /<seq_id>/residue   (L, 1024) float16  per-residue embeddings
    /<seq_id>/mean      (1024,)  float16  mean-pooled per-protein embedding

Incremental (incremental.py): an existing --out file is appended to, never rewritten. Only proteins that are not in it yet, or whose sequence changed since their vector was made,
are embedded; every group carries the md5 of its sequence (attribute input_md5). For the archived files, which have no fingerprints, give --ledger (a FASTA / CSV of the sequences
they were made from): matching vectors are kept and stamped, the others embedded again. --overwrite starts the file again. The model is not even loaded when nothing is to be done.
    python embed_prostt5.py --fasta wanted.fasta --out embeddings_SAM_ProstT5.h5 --no-per-residue --ledger archived_sequences.csv

Notes
-----
- ProstT5 is bidirectional; for encoding AA -> embedding we prepend "<AA2fold>"
  and feed uppercase sequence. Non-standard residues (U, Z, O, B) are mapped to X.
- Uses the encoder only (T5EncoderModel). If you instead want to *translate*
  AA -> 3Di tokens, use AutoModelForSeq2SeqLM.from_pretrained(...).generate(...).
"""

from __future__ import annotations

import argparse
import re
import time
from pathlib import Path

import h5py
import numpy as np
import torch
from transformers import T5EncoderModel, T5Tokenizer

import incremental as inc

MODEL_NAME = "Rostlab/ProstT5"


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    sid, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if sid is not None:
                    seqs[sid] = "".join(buf)
                sid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if sid is not None:
            seqs[sid] = "".join(buf)
    return seqs


def prepare(seq: str) -> str:
    # Map non-standard amino acids to X, space-separate residues, prepend prefix.
    seq = re.sub(r"[UZOB]", "X", seq.upper())
    return "<AA2fold> " + " ".join(seq)


def batch_by_length(items: list[tuple[str, str]], max_residues: int, max_seqs: int):
    # Sort by length so padding is minimal; then greedy-pack into batches.
    items = sorted(items, key=lambda kv: len(kv[1]))
    batch: list[tuple[str, str]] = []
    cur_max = 0
    for sid, seq in items:
        L = len(seq)
        proj = max(cur_max, L) * (len(batch) + 1)
        if batch and (proj > max_residues or len(batch) >= max_seqs):
            yield batch
            batch, cur_max = [], 0
        batch.append((sid, seq))
        cur_max = max(cur_max, L)
    if batch:
        yield batch


@torch.inference_mode()
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path, help="HDF5 output path")
    ap.add_argument("--per-residue", action="store_true", default=True,
                    help="Store per-residue embeddings (default: on)")
    ap.add_argument("--no-per-residue", dest="per_residue", action="store_false")
    ap.add_argument("--max-residues", type=int, default=4000,
                    help="Soft cap on batch-token count (batch_size * max_len).")
    ap.add_argument("--max-seqs", type=int, default=16)
    ap.add_argument("--fp16", action="store_true", default=True,
                    help="Run model in fp16 on GPU (default: on).")
    ap.add_argument("--fp32", dest="fp16", action="store_false")
    ap.add_argument("--ledger", nargs="+", help="what the unstamped vectors already in --out were made from (FASTA, or CSV / TSV with Entry and Sequence)")
    ap.add_argument("--trust-existing", action="store_true", help="keep unstamped vectors by accession alone, without a ledger")
    ap.add_argument("--overwrite", action="store_true", help="start --out again instead of appending")
    args = ap.parse_args()

    seqs = read_fasta(args.fasta)
    wanted = {sid: inc.md5_text(seq) for sid, seq in seqs.items()}
    if args.overwrite or not args.out.exists():
        pl = {"todo": {sid: "new" for sid in wanted}, "keep": [], "stamp": {}, "extra": []}
    else:
        pl = inc.plan(wanted, inc.stored_fingerprints(args.out), inc.read_ledger(args.ledger), args.trust_existing)
    inc.report("ProstT5", pl, len(seqs))
    if not pl["todo"]:
        inc.apply_stamps(args.out, pl["stamp"])
        print("[done] nothing to embed", flush=True)
        return

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[device] {device}")

    print(f"[load] {MODEL_NAME}")
    tokenizer = T5Tokenizer.from_pretrained(MODEL_NAME, do_lower_case=False)
    model = T5EncoderModel.from_pretrained(MODEL_NAME).to(device).eval()
    if args.fp16 and device.type == "cuda":
        model = model.half()
    print(f"[load] params = {sum(p.numel() for p in model.parameters())/1e6:.1f} M")

    items = [(sid, seqs[sid]) for sid in pl["todo"]]
    print(f"[data] {len(items)} sequences to embed  (avg len = {np.mean([len(s) for _, s in items]):.0f})")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    h5 = h5py.File(args.out, "w" if args.overwrite else "a")
    inc.check_config(h5, {"model": MODEL_NAME, "precision": "fp16" if (args.fp16 and device.type == "cuda") else "fp32"})

    done = 0
    t0 = time.time()
    for batch in batch_by_length(items, args.max_residues, args.max_seqs):
        ids = [b[0] for b in batch]
        raws = [b[1] for b in batch]
        prepped = [prepare(s) for s in raws]
        enc = tokenizer.batch_encode_plus(
            prepped, add_special_tokens=True, padding="longest", return_tensors="pt"
        )
        input_ids = enc["input_ids"].to(device)
        attn = enc["attention_mask"].to(device)

        out = model(input_ids=input_ids, attention_mask=attn)
        # last_hidden_state: (B, T, 1024). T includes <AA2fold> prefix + EOS.
        emb = out.last_hidden_state.float().cpu().numpy()

        for i, (sid, seq) in enumerate(zip(ids, raws)):
            L = len(seq)
            # prefix occupies token index 0; residues are tokens [1 : 1+L]
            residue_emb = emb[i, 1 : 1 + L].astype(np.float16)
            mean_emb = residue_emb.mean(axis=0).astype(np.float16)
            inc.replace_slot(h5, sid)
            grp = h5.create_group(sid)
            grp.attrs[inc.ATTR] = wanted[sid]
            if args.per_residue:
                grp.create_dataset("residue", data=residue_emb, compression="gzip")
            grp.create_dataset("mean", data=mean_emb)
            grp.attrs["length"] = L

        done += len(batch)
        print(f"[progress] {done}/{len(items)}  "
              f"({done / max(time.time() - t0, 1e-9):.1f} seq/s)")

    h5.close()
    inc.apply_stamps(args.out, pl["stamp"])
    print(f"[done] {args.out}: embedded {len(items)}, kept {len(pl['keep'])}")


if __name__ == "__main__":
    main()
