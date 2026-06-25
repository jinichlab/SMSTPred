#!/usr/bin/env python3
"""
ProstT5 self-attention vs classifier task-importance (IG), per residue.

Companion to integrated_gradients.py. IG tells us which residues drive the
*classifier's* substrate call. This asks a different question: which residues
does *ProstT5 itself* internally attend to? If the model's attention concentrates
on the same catalytic / SAM-binding residues that IG flags, that's independent
evidence the structure-aware PLM has learned those residues are special.

CAVEATS (state these in any writeup)
------------------------------------
- ProstT5 attention reflects its PRETRAINING objective (masked-span / 3Di), NOT
  our substrate task. High attention => salient for the model's own representation
  (often = conserved/structurally central residues), not a causal task explanation.
- "Attention is not explanation" -- this is corroborative, not proof.

WHAT IT COMPUTES (per showcase protein)
---------------------------------------
- attn_received : mean over all layers & heads of the attention each residue
  RECEIVES from all query residues (special tokens stripped). The simple,
  robust "what is the model focusing on" measure.
- attn_rollout  : attention rollout (Abnar & Zuidema 2020) received score, which
  accounts for multi-layer information mixing.
- overlays both against the IG track (ig_output_logit/) and marks catalytic
  motifs found by regex (SAM glycine loop, DPPY-type, PC-type catalytic).
- reports Spearman(attn_received, |IG|) and where the motif residues rank.

Env: classifier_env (GPU). Output: attn_output/.
"""
import argparse
import glob
import re
from pathlib import Path

import h5py  # noqa: F401  (kept for parity; not strictly needed)
import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr
from transformers import T5EncoderModel, T5Tokenizer

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
MODEL_NAME = "Rostlab/ProstT5"
SUBSTRATE_COLORS = {"Other": "#e6550d", "Protein": "#238b45",
                    "RNA": "#cb181d", "DNA": "#6a51a3"}

# Catalytic / cofactor sequence motifs (regex) to mark as "known important".
MOTIFS = [
    ("SAM G-loop", re.compile(r"G.G.G")),        # Rossmann glycine-rich Motif I
    ("N-cat DPPY", re.compile(r"[DNSH]PP[YFWH]")),  # N6-adenine/N4-cytosine catalytic
    ("C5-cat PC", re.compile(r"PC[QGASN]")),     # C5-cytosine catalytic Cys (PCQ-type)
]


def prepare(seq):
    seq = re.sub(r"[UZOB]", "X", seq.upper())
    return "<AA2fold> " + " ".join(seq)


def read_fasta(path):
    seqs, sid, buf = {}, None, []
    for ln in open(path):
        ln = ln.rstrip()
        if not ln:
            continue
        if ln.startswith(">"):
            if sid:
                seqs[sid] = "".join(buf)
            sid, buf = ln[1:].split()[0], []
        else:
            buf.append(ln)
    if sid:
        seqs[sid] = "".join(buf)
    return seqs


@torch.inference_mode()
def attention_profiles(model, tok, seq, device):
    """Return (attn_received, attn_rollout) per residue (length L)."""
    L = len(seq)
    enc = tok.batch_encode_plus([prepare(seq)], add_special_tokens=True,
                                padding="longest", return_tensors="pt")
    ids = enc["input_ids"].to(device)
    mask = enc["attention_mask"].to(device)
    out = model(input_ids=ids, attention_mask=mask, output_attentions=True)
    # residues are tokens [1 : 1+L]  (token 0 = <AA2fold> prefix; then EOS)
    res = slice(1, 1 + L)

    received_layers = []
    rollout = np.eye(L, dtype=np.float64)
    for layer_attn in out.attentions:                 # (1, H, T, T)
        a = layer_attn[0].mean(0).float().cpu().numpy()  # mean over heads -> (T,T)
        a = a[res, res]                                # residues only (q x k)
        received_layers.append(a.mean(axis=0))         # attention each key receives
        # rollout: add residual, renormalize rows, propagate
        a_res = a + np.eye(L)
        a_res = a_res / a_res.sum(axis=1, keepdims=True)
        rollout = a_res @ rollout
    attn_received = np.mean(received_layers, axis=0)
    attn_rollout = rollout.mean(axis=0)                # received score after mixing
    return attn_received, attn_rollout


def find_motifs(seq):
    """Return list of (name, start, end1based_inclusive) and a per-residue mask."""
    hits, mask = [], np.zeros(len(seq), dtype=bool)
    for name, rx in MOTIFS:
        for m in rx.finditer(seq):
            hits.append((name, m.start() + 1, m.end()))   # 1-based inclusive
            mask[m.start():m.end()] = True
    return hits, mask


def load_ig(entry):
    """Find the IG residue-attr CSV for this entry (logit-space run)."""
    cands = glob.glob(str(HERE / "ig_output_logit" / f"{entry}__*_residue_attr.csv"))
    if not cands:
        return None, None
    lab = Path(cands[0]).stem.split("__")[1].replace("_residue_attr", "")
    return pd.read_csv(cands[0]), lab


def pct_rank(arr, mask):
    """Mean percentile rank (0-100) of the masked positions within arr."""
    if mask.sum() == 0:
        return np.nan
    order = arr.argsort().argsort()           # 0..L-1 rank
    return float(100.0 * order[mask].mean() / (len(arr) - 1))


def plot_overlay(entry, seq, ig, attn, motif_hits, label, color, corr, out_png):
    L = len(seq)
    x = np.arange(1, L + 1)
    fig, ax1 = plt.subplots(figsize=(max(8, L / 22), 3.6))
    pos = np.clip(ig, 0, None); neg = np.clip(ig, None, 0)
    ax1.bar(x, pos, width=1.0, color=color, linewidth=0, alpha=0.85)
    ax1.bar(x, neg, width=1.0, color="#bdbdbd", linewidth=0)
    ax1.axhline(0, color="k", lw=0.5)
    ax1.set_ylabel(f"IG -> {label}", color=color)
    ax1.set_xlabel("residue position")
    ax1.set_xlim(0.5, L + 0.5)
    ax1.spines["top"].set_visible(False)

    ax2 = ax1.twinx()
    an = attn / (attn.max() + 1e-12)
    ax2.plot(x, an, color="#111111", lw=1.1, label="ProstT5 attention received")
    ax2.set_ylabel("attention received (norm)", color="#111111")
    ax2.set_ylim(0, 1.05)
    ax2.spines["top"].set_visible(False)

    for name, s, e in motif_hits:
        ax1.axvspan(s - 0.5, e + 0.5, color="#ffd92f", alpha=0.35, zorder=0)
        ax1.annotate(name, ((s + e) / 2, ax1.get_ylim()[1]),
                     ha="center", va="top", fontsize=6.5, color="#7a6a00")
    ax1.set_title(f"{entry}  {label}   Spearman(attn, |IG|) = {corr:.2f}", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", type=Path, default=HERE / "showcase.fasta")
    ap.add_argument("--out-dir", type=Path, default=HERE / "attn_output")
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[device] {device}")
    tok = T5Tokenizer.from_pretrained(MODEL_NAME, do_lower_case=False)
    model = T5EncoderModel.from_pretrained(MODEL_NAME).to(device).eval()
    if device.type == "cuda":
        model = model.half()
    print(f"[load] {MODEL_NAME}  ({sum(p.numel() for p in model.parameters())/1e6:.0f}M)")

    seqs = read_fasta(args.fasta)
    summary = []
    for entry, seq in seqs.items():
        attn_recv, attn_roll = attention_profiles(model, tok, seq, device)
        ig_df, label = load_ig(entry)
        motif_hits, motif_mask = find_motifs(seq)

        if ig_df is None:
            print(f"  [warn] no IG csv for {entry}; attention-only")
            ig = np.zeros(len(seq)); label, color = "?", "#333333"; corr = np.nan
        else:
            ig = ig_df["ig"].values
            color = SUBSTRATE_COLORS.get(label, "#333333")
            corr = spearmanr(attn_recv, np.abs(ig)).statistic

        pd.DataFrame({
            "position": np.arange(1, len(seq) + 1), "residue": list(seq),
            "ig": ig, "attn_received": attn_recv, "attn_rollout": attn_roll,
            "in_motif": motif_mask.astype(int),
        }).to_csv(args.out_dir / f"{entry}_attn_vs_ig.csv", index=False)

        plot_overlay(entry, seq, ig, attn_recv, motif_hits, label, color, corr,
                     args.out_dir / f"{entry}_attn_overlay.png")

        row = {"Entry": entry, "label": label,
               "spearman_attn_absIG": round(float(corr), 3) if corr == corr else np.nan,
               "n_motif_res": int(motif_mask.sum()),
               "motif_attn_pctile": round(pct_rank(attn_recv, motif_mask), 1),
               "motif_IG_pctile": round(pct_rank(np.abs(ig), motif_mask), 1),
               "motifs": "; ".join(f"{n}@{s}-{e}" for n, s, e in motif_hits)}
        summary.append(row)
        print(f"  {entry} ({label}): Spearman(attn,|IG|)={row['spearman_attn_absIG']}  "
              f"motif attn pct={row['motif_attn_pctile']}  IG pct={row['motif_IG_pctile']}  "
              f"[{row['motifs']}]")

    pd.DataFrame(summary).to_csv(args.out_dir / "attn_vs_ig_summary.csv", index=False)
    print(f"[done] {len(summary)} proteins -> {args.out_dir}")


if __name__ == "__main__":
    main()
