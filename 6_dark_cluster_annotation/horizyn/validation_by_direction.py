#!/usr/bin/env python3
"""Does Horizyn put a SAM-methylation reaction first for the 3,478 dark proteins? Forward direction only versus either direction.

Rhea lists every reaction twice, once per direction (master n; entries n+1 and n+2). Horizyn scores the two entries almost
identically, so asking whether the single best-scoring entry is a SAM-methylation reaction in its forward direction undercounts:
the same reaction can win in the reverse direction. This script
  1. pairs each entry with its twin by matching the two sides of the reaction SMILES (direction-collapsed reaction),
  2. shows how alike the twins are (embedding cosine, score difference),
  3. reports the top-1 hit rate both ways with the chance level under each definition.

Inputs (from the Zenodo record, in data/ or $SMSTPRED_DATA): horizyn_dark_scores.npz, horizyn_rhea_reaction_embeddings.npz;
and rhea_reaction_clusters.csv (in this folder).
Outputs: validation_by_direction.txt, validation_by_direction.csv (one row per dark protein).
"""
from pathlib import Path
import os

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = Path(os.environ.get("SMSTPRED_DATA", HERE.parents[1] / "data"))

z = np.load(DATA / "horizyn_dark_scores.npz", allow_pickle=True)
S, acc, clus, rid, is_sam = z["scores"], z["accessions"], z["clusters"], z["rhea_ids"], z["is_sam_methylation"]
emb = np.load(DATA / "horizyn_rhea_reaction_embeddings.npz", allow_pickle=True)
assert (emb["rhea_ids"] == rid).all()
E = emb["embeds"].astype(np.float64)
rc = pd.read_csv(HERE / "rhea_reaction_clusters.csv").set_index("rhea_lr_id").loc[rid]

# direction-collapsed reaction id: sorted fragments of both sides of the SMILES
side = lambda x: ".".join(sorted(x.split(".")))
keys = [tuple(sorted(side(q) for q in smi.split(">>"))) for smi in rc.reaction_smiles]
kid = {k: i for i, k in enumerate(dict.fromkeys(keys))}
rxn = np.array([kid[k] for k in keys])
n_rxn = rxn.max() + 1
sam_rxn = np.isin(rxn, np.unique(rxn[is_sam]))              # every entry (either direction) of a SAM-methylation reaction

# twins
members = pd.Series(np.arange(len(rid))).groupby(rxn).apply(list)
pairs = [m for m in members if len(m) == 2]
a, b = np.array([m[0] for m in pairs]), np.array([m[1] for m in pairs])
cos = (E[a] * E[b]).sum(1) / (np.linalg.norm(E[a], axis=1) * np.linalg.norm(E[b], axis=1))
diff = S[:, b] - S[:, a]                                    # higher-ID twin minus lower-ID twin, per protein and reaction
lines = []
w = lines.append
w(f"Rhea entries scored: {len(rid):,}; direction-collapsed reactions: {n_rxn:,} ({len(pairs):,} with exactly two entries)")
w(f"SAM-methylation entries: forward-flagged {int(is_sam.sum()):,}; with the reverse twin {int(sam_rxn.sum()):,} entries "
  f"= {len(np.unique(rxn[sam_rxn])):,} reactions")
w(f"Twin entries, embedding cosine: median {np.median(cos):.3f}, 5th percentile {np.percentile(cos, 5):.3f}")
w(f"Twin entries, score difference (higher-ID entry minus lower-ID entry) over all protein-reaction pairs: mean {diff.mean():+.4f}, "
  f"mean absolute {np.abs(diff).mean():.4f}, SD {diff.std():.4f}; higher-ID twin scores higher in {100 * (diff > 0).mean():.1f}%, "
  f"lower in {100 * (diff < 0).mean():.1f}%, exactly equal in {100 * (diff == 0).mean():.1f}% of pairs")

top1 = S.argmax(1)
hit_fwd = is_sam[top1]
hit_any = sam_rxn[top1]
n = len(top1)
chance_fwd = is_sam.mean()
chance_any = sam_rxn.mean()
w("")
w(f"Dark proteins scored: {n:,}")
w(f"Top-1 is a SAM-methylation reaction, forward direction only : {hit_fwd.sum():,}/{n:,} = {100 * hit_fwd.mean():.1f}%   (chance {100 * chance_fwd:.1f}%)")
w(f"Top-1 is a SAM-methylation reaction, either direction       : {hit_any.sum():,}/{n:,} = {100 * hit_any.mean():.1f}%   (chance {100 * chance_any:.1f}%)")
k = 100
topk = np.argpartition(-S, k, axis=1)[:, :k]
w(f"Top-100 contains at least one SAM-methylation reaction (forward only / either direction): "
  f"{100 * is_sam[topk].any(1).mean():.1f}% / {100 * sam_rxn[topk].any(1).mean():.1f}% "
  f"(chance {100 * (1 - (1 - chance_fwd) ** k):.0f}% / {100 * (1 - (1 - chance_any) ** k):.0f}%): uninformative, as the chance levels are near 100%")

# proteins whose top hit is not a SAM-methylation reaction, by cluster
miss = pd.DataFrame({"cluster": clus[~hit_any]}).value_counts().rename("n_missed").reset_index()
w("")
w(f"Proteins missed (either direction): {int((~hit_any).sum()):,}; clusters with the most misses: " +
  ", ".join(f"c{int(c)} ({int(m)})" for c, m in zip(miss.cluster[:6], miss.n_missed[:6])))

pd.DataFrame({"cluster": clus, "accession": acc, "top1_rhea_entry": rid[top1], "top1_score": S[np.arange(n), top1],
              "top1_sam_forward_only": hit_fwd, "top1_sam_either_direction": hit_any}).to_csv(HERE / "validation_by_direction.csv", index=False)
(HERE / "validation_by_direction.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
