#!/usr/bin/env python3
"""
Stage 3 / Script 2 - predict SAM substrate type for new proteins.

Loads the deployment model built by 2_training/train_deployment_ensemble.py
(6-member ProstT5 MLP ensemble + frozen per-label thresholds + applicability-
domain cutoff) and scores new proteins from their ProstT5 embeddings.

Per protein it reports, for each label:
  - score_<label>   ensemble-mean probability (0-1)
  - call_<label>    1 if score >= the frozen threshold, else 0  (multi-label)
and per protein:
  - disagreement    mean over labels of the 6 members' score spread (SD); high = shaky
  - ad_maxsim       max cosine similarity to the training set (ProstT5 space)
  - in_domain       1 if ad_maxsim >= AD cutoff else 0  (0 = extrapolating)

Input: one ProstT5 embedding H5 file, or a directory of them (Entry -> /mean or
flat vector), as produced by embed_prostt5.py on the target sequences. Writes one
predictions CSV per input file to --out-dir.

Env: smstpred (`make env`) (GPU).
"""
import os
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
import torch
import torch.nn as nn

REPO = Path(__file__).resolve().parent.parent
MODEL_DIR = REPO / "2_training" / "deployment_model"
FILES_DIR = Path(os.environ.get("SMSTPRED_DATA", REPO / "data"))
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
BATCH = 4096


class MultilabelNN(nn.Module):
    def __init__(self, in_dim, out_dim, hidden):
        super().__init__()
        self.fc1 = nn.Linear(in_dim, hidden[0])
        self.fc2 = nn.Linear(hidden[0], hidden[1])
        self.fc3 = nn.Linear(hidden[1], out_dim)
        self.relu = nn.ReLU()

    def forward(self, x):
        return self.fc3(self.relu(self.fc2(self.relu(self.fc1(x)))))


def load_embeddings(path):
    """Entry -> vector, accepting either nested {Entry}/mean or a flat dataset."""
    out = {}
    with h5py.File(path, "r") as f:
        for k in f.keys():
            item = f[k]
            out[k] = np.array(item["mean"] if isinstance(item, h5py.Group) else item,
                              dtype=np.float32)
    return out


def load_model():
    cfg = json.load(open(MODEL_DIR / "model_config.json"))
    thr_j = json.load(open(MODEL_DIR / "thresholds.json"))
    labels = cfg["labels"]
    thresholds = np.array([thr_j["thresholds"][l] for l in labels], dtype=np.float32)
    models = []
    for s in range(cfg["n_ensemble"]):
        m = MultilabelNN(cfg["input_dim"], len(labels), cfg["hidden"]).to(DEVICE)
        m.load_state_dict(torch.load(MODEL_DIR / f"ensemble_member_{s}.pt", map_location=DEVICE))
        m.eval()
        models.append(m)
    return cfg, labels, thresholds, models


def load_ad_reference(cfg):
    ref = load_embeddings(FILES_DIR / cfg["ad_reference_h5"])
    X = np.vstack(list(ref.values())).astype(np.float32)
    Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-8)
    return torch.tensor(Xn, device=DEVICE)


@torch.no_grad()
def predict_file(path, labels, thresholds, models, ref, ad_cutoff, out_csv):
    emb = load_embeddings(path)
    entries = list(emb.keys())
    if not entries:
        print(f"  {Path(path).name}: empty, skipped")
        return 0
    X = np.vstack([emb[e] for e in entries]).astype(np.float32)
    thr_t = torch.tensor(thresholds, device=DEVICE)
    rows = []
    for i in range(0, len(entries), BATCH):
        xb = torch.tensor(X[i:i + BATCH], dtype=torch.float32, device=DEVICE)
        member_probs = torch.stack([torch.sigmoid(m(xb)) for m in models])  # (M,B,L)
        mean = member_probs.mean(0)                       # (B,L)
        disagree = member_probs.std(0).mean(1)            # (B,)
        calls = (mean >= thr_t).int()
        xn = xb / (xb.norm(dim=1, keepdim=True) + 1e-8)
        maxsim = (xn @ ref.T).max(1).values               # (B,)
        in_dom = (maxsim >= ad_cutoff).int()
        mean, calls = mean.cpu().numpy(), calls.cpu().numpy()
        disagree, maxsim, in_dom = (disagree.cpu().numpy(), maxsim.cpu().numpy(),
                                    in_dom.cpu().numpy())
        for j, e in enumerate(entries[i:i + BATCH]):
            r = {"Entry": e}
            r.update({f"score_{lab}": round(float(mean[j, li]), 4)
                      for li, lab in enumerate(labels)})
            r.update({f"call_{lab}": int(calls[j, li]) for li, lab in enumerate(labels)})
            r["disagreement"] = round(float(disagree[j]), 4)
            r["ad_maxsim"] = round(float(maxsim[j]), 4)
            r["in_domain"] = int(in_dom[j])
            rows.append(r)
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    return len(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="ProstT5 embedding H5 file or directory")
    ap.add_argument("--pattern", default="*.h5", help="glob when --input is a directory")
    ap.add_argument("--out-dir", default="predictions_output")
    args = ap.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg, labels, thresholds, models = load_model()
    print(f"device={DEVICE} | labels={labels}")
    print(f"thresholds={ {l: round(float(t), 3) for l, t in zip(labels, thresholds)} }")
    print(f"AD reference={cfg['ad_reference_h5']}  cutoff={cfg['ad_cutoff']:.4f}")
    ref = load_ad_reference(cfg)
    ad_cutoff = float(cfg["ad_cutoff"])

    inp = Path(args.input)
    files = [inp] if inp.is_file() else sorted(inp.glob(args.pattern))
    print(f"{len(files)} prediction file(s)")
    total = 0
    for path in files:
        out_csv = out_dir / f"{path.stem}_predictions.csv"
        total += predict_file(path, labels, thresholds, models, ref, ad_cutoff, out_csv)
        print(f"  {path.name} -> {out_csv.name}")
    print(f"done: {total} proteins across {len(files)} file(s)")


if __name__ == "__main__":
    main()
