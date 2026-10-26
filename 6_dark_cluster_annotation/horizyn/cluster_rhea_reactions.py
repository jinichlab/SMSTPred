#!/usr/bin/env python3
"""Cluster all 36,552 Rhea reactions in HORIZYN'S OWN query-embedding space
(512-dim, L2-normalized -- the space the model was trained so that a reaction
sits near the enzymes that catalyze it), using HDBSCAN, matching the
clustering methodology already used elsewhere in this project for proteins.

Rationale: c758's dominant vote looked "heterogeneous" (H=0.66) purely because
Rhea encodes each isoprenoid-chain-length variant of menaquinone as a SEPARATE
reaction ID, fragmenting one real chemistry into ~10 IDs. Clustering reactions
by embedding similarity (not raw Rhea ID) should merge exactly that kind of
near-duplicate split, and may reveal other real orphan-activity clusters that
the raw-ID analysis missed or fragmented, as well as give cleaner "family"
labels and orphan-ness scores.

Outputs:
  rhea_reaction_clusters.csv   rhea_lr_id, reaction_cluster, is_sam_methylation,
                                has_enzyme, equation (where available)
  rhea_reaction_embeds.npz     the (n, 512) embeddings + ids, for reuse
"""
import os
import sys
from os.path import join, dirname, abspath

import numpy as np
import pandas as pd
import torch
import umap
from sklearn.cluster import HDBSCAN

# Matches the UMAP convention already established in this project for embedding
# clustering (repo/5_umap/global_reumap_prostt5.py): cosine metric, n_neighbors=15,
# min_dist=0.1, random_state=42, 2D. HDBSCAN min_cluster_size/min_samples are NOT
# reused from that script (50/10, tuned for large enzyme families among 112,000
# proteins) -- here we want small chemical "families" (a handful of chain-length
# variants), so a much smaller scale is used instead.
RANDOM_SEED, N_NEIGHBORS, MIN_DIST, UMAP_METRIC = 42, 15, 0.1, "cosine"
HDBSCAN_MIN_CLUSTER_SIZE, HDBSCAN_MIN_SAMPLES = 5, 3

HERE = dirname(abspath(__file__))
sys.path.insert(0, HERE)

from horizyn.config import load_config  # noqa: E402
from horizyn.datasets.base import BaseDataset  # noqa: E402
from horizyn.datasets.collection import MergeDataset  # noqa: E402
from horizyn.datasets.fingerprints import DRFPFingerprintDataset, RDKitPlusFingerprintDataset  # noqa: E402
from horizyn.datasets.transform import ConcatTensorTransform  # noqa: E402
from horizyn.lightning_module import HorizynLitModule  # noqa: E402

CHECKPOINT = join(HERE, "checkpoints", "horizyn_v1_0_inf.ckpt")
CONFIG = join(HERE, "configs", "sota.yaml")
INPUTS = os.environ.get("SMSTPRED_HORIZYN_INPUTS", join(HERE, "inputs"))     # rhea-reaction-smiles.tsv and sam_methylation_reactions_smiles.csv
RHEA_MASTER_TSV = join(INPUTS, "rhea-reaction-smiles.tsv")
SAM_METHYLATION_CSV = join(INPUTS, "sam_methylation_reactions_smiles.csv")

OUT_CSV = join(HERE, "rhea_reaction_clusters.csv")
OUT_NPZ = join(HERE, "rhea_reaction_embeds.npz")

FP_LOG_EVERY = 2000


def fingerprint_all_reactions(config):
    rhea = pd.read_csv(RHEA_MASTER_TSV, sep="\t", header=None, names=["id", "smiles"])
    sam = pd.read_csv(SAM_METHYLATION_CSV)
    sam_lr_ids = set(sam.lr_id.astype(int))
    has_enzyme_by_lr = dict(zip(sam.lr_id.astype(int), sam.has_enzyme))

    keys_all = [f"rxn_{i}" for i in range(len(rhea))]
    data_all = [{"reaction_smiles": s} for s in rhea.smiles]
    reactions = BaseDataset(keys=keys_all, array_data=data_all)

    rdkit_fp = RDKitPlusFingerprintDataset(
        reaction_dataset=reactions, vec_dim=config.data.get("rdkit_fp_dim", 1024),
        mol_fp_type="morgan", rxn_fp_type="struct", use_chirality=True,
        standardize=config.data.get("standardize_reactions", True),
        standardize_hypervalent=config.data.get("standardize_hypervalent", True),
        standardize_remove_hs=config.data.get("standardize_remove_hs", True),
        standardize_kekulize=config.data.get("standardize_kekulize", False),
        standardize_uncharge=config.data.get("standardize_uncharge", True),
        standardize_metals=config.data.get("standardize_metals", True),
    )
    drfp_fp = DRFPFingerprintDataset(
        reaction_dataset=reactions, vec_dim=config.data.get("drfp_dim", 1024),
        radius=3, rings=True,
        standardize=config.data.get("standardize_reactions", True),
        standardize_hypervalent=config.data.get("standardize_hypervalent", True),
        standardize_remove_hs=config.data.get("standardize_remove_hs", True),
        standardize_kekulize=config.data.get("standardize_kekulize", False),
        standardize_uncharge=config.data.get("standardize_uncharge", True),
        standardize_metals=config.data.get("standardize_metals", True),
    )
    merged = MergeDataset(datasets={"rdkit": rdkit_fp, "drfp": drfp_fp}, add_prefix=False)
    merged.append_transforms(ConcatTensorTransform(labels=["rdkit", "drfp"], dim=0))

    fps, ids, is_sam_l, has_enz_l, smiles_l = [], [], [], [], []
    n_fail = 0
    for i, (k, rhea_id, smi) in enumerate(zip(keys_all, rhea.id, rhea.smiles)):
        try:
            fps.append(merged[k])
            ids.append(int(rhea_id))
            is_sam_l.append(int(rhea_id) in sam_lr_ids)
            has_enz_l.append(bool(has_enzyme_by_lr.get(int(rhea_id), False)))
            smiles_l.append(smi)
        except Exception:
            n_fail += 1
        if (i + 1) % FP_LOG_EVERY == 0:
            print(f"  fingerprinted {i+1}/{len(rhea)}  (failed so far: {n_fail})", flush=True)
    print(f"[cluster] {len(fps)} fingerprinted OK, {n_fail} skipped", flush=True)
    return torch.stack(fps), np.array(ids), np.array(is_sam_l), np.array(has_enz_l), smiles_l


def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[cluster] device={device}", flush=True)
    config = load_config(CONFIG)
    model = HorizynLitModule.load_from_checkpoint(CHECKPOINT, map_location=device)
    model.eval().to(device)

    print("[cluster] fingerprinting all Rhea reactions...", flush=True)
    fps, rhea_ids, is_sam, has_enzyme, smiles = fingerprint_all_reactions(config)
    fps = fps.to(device)
    with torch.no_grad():
        embeds = model.model.query_encoder(fps.float()).cpu().numpy()
    print(f"[cluster] embeddings shape={embeds.shape}", flush=True)

    print("[cluster] running UMAP (2D, cosine, matching project convention)...", flush=True)
    coords = umap.UMAP(n_components=2, n_neighbors=N_NEIGHBORS, min_dist=MIN_DIST,
                        metric=UMAP_METRIC, random_state=RANDOM_SEED, verbose=True).fit_transform(embeds)
    print("[cluster] running HDBSCAN on UMAP coords...", flush=True)
    clusterer = HDBSCAN(min_cluster_size=HDBSCAN_MIN_CLUSTER_SIZE, min_samples=HDBSCAN_MIN_SAMPLES)
    labels = clusterer.fit_predict(coords)
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = int((labels == -1).sum())
    print(f"[cluster] {n_clusters} reaction clusters, {n_noise} noise points "
          f"({n_noise/len(labels)*100:.1f}%)", flush=True)

    np.savez_compressed(OUT_NPZ, embeds=embeds, umap_coords=coords, rhea_ids=rhea_ids, labels=labels)
    df = pd.DataFrame({
        "rhea_lr_id": rhea_ids, "reaction_cluster": labels,
        "is_sam_methylation": is_sam, "has_enzyme": has_enzyme, "reaction_smiles": smiles,
    })
    df.to_csv(OUT_CSV, index=False)
    print(f"[cluster] wrote {OUT_CSV} and {OUT_NPZ}", flush=True)
    print("[cluster] DONE", flush=True)


if __name__ == "__main__":
    main()
