#!/usr/bin/env python3
"""SMSTPred ready classifier: protein sequences in (FASTA), substrate type of the SAM-methyltransferase out (CSV).

  python classify.py my_proteins.fasta                       -> my_proteins_predictions.csv
  python classify.py my_proteins.fasta -o results.csv --device cpu

What it does, in order: reads and checks the FASTA (sequences over 2,000 aa, empty records and non-letters are reported, never silently dropped); embeds every protein with ProtT5-XL
(Rostlab/prot_t5_xl_half_uniref50-enc, mean over the residues; downloaded from Hugging Face on the first run, about 1.2 GB); scores the embeddings with the newest trained ensemble
(get_model.py: model/LATEST.json) and writes one row per input protein. See README.md for the input format and every output column.
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from get_model import resolve  # noqa: E402

MAX_LEN = 2000                    # the longest length the model was trained and validated for
MIN_LEN = 30                      # below this a sequence is a fragment or a peptide, not a methyltransferase domain protein
PROT_T5 = "Rostlab/prot_t5_xl_half_uniref50-enc"


def read_fasta(path):
    """[(id, sequence)] in file order; ids are the first word of the header, made unique if repeated."""
    recs, sid, buf, seen = [], None, [], {}
    def flush():
        if sid is not None:
            recs.append((sid, "".join(buf)))
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            flush()
            sid = line[1:].split()[0] if line[1:].split() else "unnamed"
            seen[sid] = seen.get(sid, 0) + 1
            if seen[sid] > 1:
                sid = f"{sid}_{seen[sid]}"
            buf = []
        elif line and sid is not None:
            buf.append("".join(line.split()).upper().replace("-", "").replace("*", ""))
    flush()
    return recs


class MultilabelNN(nn.Module):
    def __init__(self, in_dim, out_dim, hidden):
        super().__init__()
        self.fc1, self.fc2, self.fc3 = nn.Linear(in_dim, hidden[0]), nn.Linear(hidden[0], hidden[1]), nn.Linear(hidden[1], out_dim)
        self.relu = nn.ReLU()

    def forward(self, x):
        return self.fc3(self.relu(self.fc2(self.relu(self.fc1(x)))))


@torch.inference_mode()
def embed(seqs, device, max_residues=4000, max_seqs=64):
    """{id: 1024-d vector}: ProtT5-XL encoder, mean over the residues (the pipeline's training embedding)."""
    from transformers import T5EncoderModel, T5Tokenizer
    half = device.type == "cuda"
    tok = T5Tokenizer.from_pretrained(PROT_T5, do_lower_case=False)
    model = T5EncoderModel.from_pretrained(PROT_T5, torch_dtype=torch.float16 if half else torch.float32).to(device).eval()
    items = sorted(seqs.items(), key=lambda kv: -len(kv[1]))
    out, batch, longest = {}, [], 0
    def run(b):
        enc = tok([" ".join(re.sub(r"[UZOB]", "X", s)) for _, s in b], add_special_tokens=True, padding="longest", return_tensors="pt").to(device)
        h = model(input_ids=enc.input_ids, attention_mask=enc.attention_mask).last_hidden_state
        for i, (sid, s) in enumerate(b):
            out[sid] = h[i, :len(s)].mean(dim=0).float().cpu().numpy()
        print(f"  embedded {len(out)}/{len(items)}", flush=True)
    for sid, s in items:
        if batch and (max(longest, len(s)) * (len(batch) + 1) > max_residues or len(batch) >= max_seqs):
            run(batch)
            batch, longest = [], 0
        batch.append((sid, s))
        longest = max(longest, len(s))
    if batch:
        run(batch)
    return out


@torch.inference_mode()
def score(X, folder, device):
    """The ensemble of the model folder on an (n, 1024) matrix: mean probability, calls, disagreement, applicability domain."""
    import h5py
    cfg = json.load(open(folder / "model_config.json"))
    labels = cfg["labels"]
    thr = torch.tensor([json.load(open(folder / "thresholds.json"))["thresholds"][l] for l in labels], device=device)
    models = []
    for s in range(cfg["n_ensemble"]):
        m = MultilabelNN(cfg["input_dim"], len(labels), cfg["hidden"]).to(device)
        m.load_state_dict(torch.load(folder / f"ensemble_member_{s}.pt", map_location=device))
        models.append(m.eval())
    with h5py.File(folder / cfg["ad_reference_path"] if cfg.get("ad_reference_path") else folder / cfg["ad_reference_h5"], "r") as f:
        ref = np.vstack([np.array(f[k]["mean"] if hasattr(f[k], "keys") else f[k], dtype=np.float32) for k in f.keys()])
    ref = torch.tensor(ref / (np.linalg.norm(ref, axis=1, keepdims=True) + 1e-8), device=device)
    xb = torch.tensor(X, dtype=torch.float32, device=device)
    probs = torch.stack([torch.sigmoid(m(xb)) for m in models])
    mean = probs.mean(0)
    xn = xb / (xb.norm(dim=1, keepdim=True) + 1e-8)
    maxsim = (xn @ ref.T).max(1).values
    return (labels, mean.cpu().numpy(), (mean >= thr).int().cpu().numpy(), probs.std(0).mean(1).cpu().numpy(), maxsim.cpu().numpy(), float(cfg["ad_cutoff"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", type=Path)
    ap.add_argument("-o", "--out", type=Path)
    ap.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    ap.add_argument("--model-dir", help="a trained ensemble folder (default: the newest model of this repository, model/LATEST.json)")
    ap.add_argument("--save-embeddings", type=Path, help="also keep the ProtT5 vectors as an .npz")
    a = ap.parse_args()
    out = a.out or a.fasta.with_name(a.fasta.stem + "_predictions.csv")
    device = torch.device("cuda" if (a.device == "auto" and torch.cuda.is_available()) or a.device == "cuda" else "cpu")
    folder, info = resolve(a.model_dir)
    print(f"model: release {info.get('release')} ({folder}); device {device}", flush=True)

    recs = read_fasta(a.fasta)
    if not recs:
        sys.exit(f"no sequences found in {a.fasta}: is it FASTA (a '>' header line, then the sequence)?")
    rows, todo = [], {}
    for sid, s in recs:
        bad = sorted(set(re.sub(r"[ACDEFGHIKLMNPQRSTVWYUZOBX]", "", s)))
        if not s:
            status = "empty sequence"
        elif len(s) < MIN_LEN:
            status = f"not scored: {len(s)} aa is shorter than {MIN_LEN} aa"
        elif len(s) > MAX_LEN:
            status = f"not scored: {len(s)} aa is longer than {MAX_LEN} aa"
        elif bad:
            status = f"not scored: unexpected characters {''.join(bad)}"
        else:
            status = "scored"
            todo[sid] = s
        rows.append({"Entry": sid, "length": len(s), "status": status})
    print(f"{len(recs)} sequences read; {len(todo)} can be scored; {len(recs) - len(todo)} cannot (see the status column)", flush=True)
    res = pd.DataFrame(rows)
    if todo:
        vecs = embed(todo, device)
        ids = list(todo)
        labels, mean, calls, dis, sim, cut = score(np.vstack([vecs[i] for i in ids]), folder, device)
        if a.save_embeddings:
            np.savez_compressed(a.save_embeddings, ids=np.array(ids), vectors=np.vstack([vecs[i] for i in ids]))
        sc = pd.DataFrame({"Entry": ids, **{f"score_{l}": np.round(mean[:, k], 4) for k, l in enumerate(labels)}, **{f"call_{l}": calls[:, k] for k, l in enumerate(labels)},
                           "disagreement": np.round(dis, 4), "ad_maxsim": np.round(sim, 4), "in_domain": (sim >= cut).astype(int)})
        sc["predicted_substrate"] = ["+".join(l for k, l in enumerate(labels) if calls[i, k]) or "none" for i in range(len(ids))]
        res = res.merge(sc, on="Entry", how="left")
    res.to_csv(out, index=False)
    n = res[res.status == "scored"]
    if len(n):
        print("\npredicted substrate type:\n" + n.predicted_substrate.value_counts().to_string())
        print(f"\ninside the applicability domain: {int(n.in_domain.sum())} of {len(n)} (the others are far from every training protein: treat their calls with caution)")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
