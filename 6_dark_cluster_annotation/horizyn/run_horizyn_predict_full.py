#!/usr/bin/env python3
"""Score ALL 3,478 dark-cluster proteins against ALL 36,552 directional Rhea
reactions (the full master reaction-SMILES file, every EC class -- not just
SAM-methylation) using the pretrained Horizyn `inf` checkpoint. No training,
no filtering to "clean small molecule": if a reaction SMILES fails to
fingerprint it is skipped and logged, everything else is scored.

The point is a VALIDATION check as much as a discovery screen: for each
protein, where do the 858 known SAM-methylation reactions (identified by
Rhea's own `lr_id`) land in that protein's full 36,552-reaction ranking? If
Horizyn has learned real chemistry, SAM-methylation reactions should rank far
above the ~2.3% baseline you'd expect by chance for these SAM-MTase proteins.

Outputs (this dir):
  horizyn_full_topk.csv       top-50 reactions per protein (compact, inspectable)
  horizyn_full_validation.csv per-protein: rank/percentile of best SAM-methylation hit
  horizyn_full_scores.npz     the raw dense (n_proteins x n_valid_reactions) score
                              matrix + protein/reaction id arrays, for later slicing
                              without re-running embeddings
"""
import os
import sys
from os.path import join, dirname, abspath

import numpy as np
import pandas as pd
import torch

HERE = dirname(abspath(__file__))
sys.path.insert(0, HERE)

from horizyn.config import load_config  # noqa: E402
from horizyn.datasets.base import BaseDataset  # noqa: E402
from horizyn.datasets.collection import MergeDataset  # noqa: E402
from horizyn.datasets.fingerprints import DRFPFingerprintDataset, RDKitPlusFingerprintDataset  # noqa: E402
from horizyn.datasets.transform import ConcatTensorTransform  # noqa: E402
from horizyn.lightning_module import HorizynLitModule  # noqa: E402

from embed_prott5 import load_prott5, embed_sequences  # noqa: E402

CHECKPOINT = join(HERE, "checkpoints", "horizyn_v1_0_inf.ckpt")
CONFIG = join(HERE, "configs", "sota.yaml")

# STAGE6 = the 6_dark_cluster_annotation/ folder of the repository
STAGE6 = os.environ.get("SMSTPRED_STAGE6", dirname(HERE))
DARK33_FASTA = join(STAGE6, "clean", "dark33.fasta")
DARK33_MEMBERS = join(STAGE6, "dark33_members.csv")
INPUTS = os.environ.get("SMSTPRED_HORIZYN_INPUTS", join(HERE, "inputs"))     # rhea-reaction-smiles.tsv and sam_methylation_reactions_smiles.csv
RHEA_MASTER_TSV = join(INPUTS, "rhea-reaction-smiles.tsv")
SAM_METHYLATION_CSV = join(INPUTS, "sam_methylation_reactions_smiles.csv")

OUT_TOPK = join(HERE, "horizyn_full_topk.csv")
OUT_VALIDATION = join(HERE, "horizyn_full_validation.csv")
OUT_NPZ = join(HERE, "horizyn_full_scores.npz")

TOP_K = 50
BATCH_SIZE_PROT = 16
FP_LOG_EVERY = 2000


def load_all_proteins():
    seqs = {}
    acc = None
    for line in open(DARK33_FASTA):
        line = line.rstrip("\n")
        if line.startswith(">"):
            acc = line[1:].split()[0]
            seqs[acc] = ""
        elif acc:
            seqs[acc] += line
    members = pd.read_csv(DARK33_MEMBERS)  # Entry, Organism, cluster
    members = members.drop_duplicates("Entry")
    cluster_of = dict(zip(members.Entry.astype(str), members.cluster))
    accs = [a for a in seqs if a in cluster_of]
    missing = [a for a in seqs if a not in cluster_of]
    if missing:
        print(f"[full] WARNING: {len(missing)} fasta accessions have no cluster in "
              f"{DARK33_MEMBERS}, dropping them", flush=True)
    return accs, [seqs[a] for a in accs], [cluster_of[a] for a in accs]


def load_all_reactions_with_fingerprints(config):
    rhea = pd.read_csv(RHEA_MASTER_TSV, sep="\t", header=None, names=["id", "smiles"])
    sam = pd.read_csv(SAM_METHYLATION_CSV)  # rhea, eq, has_enzyme, lr_id, reaction_smiles
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

    fps, valid_ids, valid_is_sam, valid_has_enz = [], [], [], []
    n_fail = 0
    for i, (k, rhea_id) in enumerate(zip(keys_all, rhea.id)):
        try:
            fps.append(merged[k])
            valid_ids.append(int(rhea_id))
            valid_is_sam.append(int(rhea_id) in sam_lr_ids)
            valid_has_enz.append(bool(has_enzyme_by_lr.get(int(rhea_id), False)))
        except Exception:
            n_fail += 1
        if (i + 1) % FP_LOG_EVERY == 0:
            print(f"  fingerprinted {i+1}/{len(rhea)}  (failed so far: {n_fail})", flush=True)

    print(f"[full] reactions: {len(rhea)} total, {len(fps)} fingerprinted OK, {n_fail} skipped "
          f"(unparseable)", flush=True)
    print(f"[full] known SAM-methylation reactions retained: {sum(valid_is_sam)} / "
          f"{len(sam_lr_ids)}", flush=True)
    return torch.stack(fps), np.array(valid_ids), np.array(valid_is_sam), np.array(valid_has_enz)


def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[full] device={device}", flush=True)

    config = load_config(CONFIG)
    print(f"[full] loading checkpoint {CHECKPOINT}", flush=True)
    model = HorizynLitModule.load_from_checkpoint(CHECKPOINT, map_location=device)
    model.eval().to(device)

    # --- proteins: all 3,478 ---
    accs, seqs, clusters = load_all_proteins()
    print(f"[full] {len(accs)} proteins to embed", flush=True)
    tokenizer, t5_model = load_prott5(device)
    prot_embeds = embed_sequences(seqs, tokenizer, t5_model, device, batch_size=BATCH_SIZE_PROT).to(device)
    del t5_model
    if device.type == "cuda":
        torch.cuda.empty_cache()
    with torch.no_grad():
        target_embeds = model.model.target_encoder(prot_embeds.float())
    print(f"[full] target_embeds shape={tuple(target_embeds.shape)}", flush=True)

    # --- reactions: the full Rhea master file, skipping unparseable entries ---
    print("[full] fingerprinting all Rhea reactions (skipping any that fail)...", flush=True)
    query_fps, rhea_ids, is_sam_meth, has_enzyme = load_all_reactions_with_fingerprints(config)
    query_fps = query_fps.to(device)
    with torch.no_grad():
        query_embeds = model.model.query_encoder(query_fps.float())
    print(f"[full] query_embeds shape={tuple(query_embeds.shape)}", flush=True)

    # --- score: cosine similarity (both sides L2-normalized) ---
    with torch.no_grad():
        scores = torch.matmul(target_embeds, query_embeds.T).cpu().numpy()  # (n_prot, n_valid_rxn)
    print(f"[full] score matrix shape={scores.shape}", flush=True)

    np.savez_compressed(OUT_NPZ, scores=scores, accessions=np.array(accs),
                         clusters=np.array(clusters), rhea_ids=rhea_ids,
                         is_sam_methylation=is_sam_meth, has_enzyme=has_enzyme)
    print(f"[full] wrote raw score matrix -> {OUT_NPZ}", flush=True)

    # --- top-K per protein ---
    topk_rows = []
    for i, (acc, cluster) in enumerate(zip(accs, clusters)):
        row_scores = scores[i]
        top_idx = np.argsort(-row_scores)[:TOP_K]
        for rank, j in enumerate(top_idx, start=1):
            topk_rows.append({
                "cluster": cluster, "accession": acc, "rank": rank,
                "rhea_lr_id": int(rhea_ids[j]), "score": float(row_scores[j]),
                "is_sam_methylation": bool(is_sam_meth[j]), "has_enzyme": bool(has_enzyme[j]),
            })
    pd.DataFrame(topk_rows).to_csv(OUT_TOPK, index=False)
    print(f"[full] wrote top-{TOP_K} per protein -> {OUT_TOPK}", flush=True)

    # --- validation: where do known SAM-methylation reactions rank, per protein? ---
    n_total = scores.shape[1]
    n_sam = int(is_sam_meth.sum())
    baseline_pct = n_sam / n_total
    val_rows = []
    for i, (acc, cluster) in enumerate(zip(accs, clusters)):
        row_scores = scores[i]
        ranking = np.argsort(-row_scores)  # best score first
        rank_of = {j: r + 1 for r, j in enumerate(ranking)}  # 1-indexed rank
        sam_positions = [rank_of[j] for j in np.where(is_sam_meth)[0]]
        best_rank = min(sam_positions)
        median_rank = float(np.median(sam_positions))
        n_in_top100 = sum(1 for r in sam_positions if r <= 100)
        val_rows.append({
            "cluster": cluster, "accession": acc,
            "best_sam_methylation_rank": best_rank,
            "best_sam_methylation_percentile": best_rank / n_total,
            "median_sam_methylation_rank": median_rank,
            "n_sam_methylation_in_top100": n_in_top100,
        })
    val_df = pd.DataFrame(val_rows)
    val_df.to_csv(OUT_VALIDATION, index=False)
    print(f"[full] wrote per-protein validation -> {OUT_VALIDATION}", flush=True)

    print(f"\n[full] === VALIDATION SUMMARY ===", flush=True)
    print(f"  total reactions scored: {n_total}  |  known SAM-methylation among them: {n_sam} "
          f"({baseline_pct*100:.2f}% -- random-chance baseline)", flush=True)
    print(f"  mean best-SAM-methylation percentile across proteins: "
          f"{val_df.best_sam_methylation_percentile.mean()*100:.2f}%  "
          f"(lower = better than chance; chance = {baseline_pct*100:.2f}%)", flush=True)
    print(f"  median best-SAM-methylation rank: {val_df.best_sam_methylation_rank.median():.0f} "
          f"of {n_total}", flush=True)
    print(f"  proteins with >=1 SAM-methylation reaction in their own top-100: "
          f"{(val_df.n_sam_methylation_in_top100 > 0).sum()} / {len(val_df)}", flush=True)
    print("[full] DONE", flush=True)


if __name__ == "__main__":
    main()
