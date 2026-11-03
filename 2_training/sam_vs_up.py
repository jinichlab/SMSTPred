#!/usr/bin/env python3
"""
Stage 2 -- SAM-trained vs UP-trained ProtT5 MLP, both scored on the SAME held-out SAM test clusters: ONE split per run (16 runs = the 16 splits of the bake-off).

  python sam_vs_up.py --iteration 3            # <release>/sam_vs_up/iterations/iter03_*
  python sam_vs_up.py --iteration 1 --smoke    # tiny end-to-end test (20 epochs, small subsets)
  sbatch run_sam_vs_up.sh (16-task array, after run_bakeoff.sh) ; then sbatch run_sam_vs_up_report.sh

The question: does training on the broad UP enzyme set (200,000 proteins) predict the substrate type of held-out SAM methyltransferases as well as training on SAM itself?
Everything but the training data is the same: ProtT5 embeddings, the same MLP, plain BCE, 2,000 epochs, early stopping (patience 80) on a validation set of unseen clusters that keeps the
best weights, a fixed 0.5 threshold, the same test proteins (the 20 % held-out SAM clusters of split i). The arms:
  ProtT5_SAM       the bake-off's MLP_ProtT5 of the same split (SAM train clusters; SAM validation clusters). Read from <release>/benchmark/iterations/iterNN_scores.npz, not retrained.
  ProtT5_UP        trained here on the UP train clusters, early-stopped on the UP validation clusters (1_clustering/up_splits.py: every UP cluster that contains a SAM test protein is
                   removed from the UP set, 15 % of the remaining UP clusters are the validation set, the rest train).
  ProtT5_SAMclean  (sensitivity) the SAM arm retrained WITHOUT the SAM train / validation proteins that share a joint SAM+UP cluster with a test protein, i.e. with the homologue exclusion
                   the UP arm has. SAM's own clusters are made on SAM alone, so some train proteins are joint-cluster mates of test proteins; the UP arm loses them, the SAM arm keeps them.
Inputs: the release folder (labels_UP.csv, labels_SAM.csv, embeddings/embeddings_UP_ProtT5.h5 and embeddings_SAM_ProtT5.h5, clusters/*). Settings: ../refresh_config.json.
Writes iterNN_detailed.csv, iterNN_per_label.csv (classifier = ProtT5_SAM / ProtT5_UP / ProtT5_SAMclean), iterNN_curves.csv.gz and iterNN_training.csv (the UP and SAMclean MLPs),
iterNN_scores.npz (their test probabilities). Resumable (--force redoes). Needs a GPU for the full run (the UP train set is about 170,000 proteins).
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

import bakeoff as b

LABELS = b.LABELS
ARMS = [a.strip() for a in os.environ.get("SAM_VS_UP_ARMS", "SAM,UP,SAMclean").split(",")]


def load_vectors(path, entries):
    """{Entry: float32 vector} for the entries that are in the H5 file (flat layout)."""
    out = {}
    with h5py.File(path, "r") as f:
        for e in entries:
            if e in f:
                out[e] = np.asarray(f[e][:], dtype=np.float32)
    return out


def stack(vecs, entries):
    ok = [e for e in entries if e in vecs]
    return ok, (np.vstack([vecs[e] for e in ok]).astype(np.float32) if ok else np.zeros((0, 1024), np.float32))


def metrics_rows(name, y_te, p_te, it, n_train, n_val, n_test, overall, per_label):
    overall[name] = {"n_train": n_train, "n_val": n_val, "n_test": n_test, **b.all_metrics(y_te, p_te, b.THRESHOLD)}
    per_label[name] = b.per_label_metrics(y_te, p_te, b.THRESHOLD)


def run_iteration(it, smoke, force):
    rel = b.release_dir()
    out_root = Path(os.environ.get("SAM_VS_UP_OUT", rel / "sam_vs_up")) / ("_smoke" if smoke else "")
    out_dir = out_root / "iterations"
    out_dir.mkdir(parents=True, exist_ok=True)
    f_det, f_pl = out_dir / f"iter{it:02d}_detailed.csv", out_dir / f"iter{it:02d}_per_label.csv"
    if f_det.exists() and not force:
        print(f"iteration {it}: {f_det.name} exists, skipping (use --force to redo)")
        return
    sp_cfg = b.CFG["splits"]
    tag = f"t{round(sp_cfg['test_frac'] * 100)}_v{round(sp_cfg['val_frac'] * 100)}"
    sam_sp = json.loads((rel / "clusters" / f"splits_tvt_sam_benchmark_id30_{tag}.json").read_text())["iterations"][it - 1]
    up_file = rel / "clusters" / f"splits_tvt_up_id30_{tag}.json"
    if not up_file.exists():
        sys.exit(f"{up_file} not found: run 1_clustering/up_splits.py (make up-splits)")
    up_sp = json.loads(up_file.read_text())["iterations"][it - 1]
    bake = Path(os.environ.get("BAKEOFF_OUT", rel / "benchmark")) / "iterations" / f"iter{it:02d}_scores.npz"
    if "SAM" in ARMS and not bake.exists():
        sys.exit(f"{bake} not found: run the bake-off (run_bakeoff.sh) first; the SAM arm is its MLP_ProtT5")
    epochs, patience = b.EPOCHS, b.PATIENCE

    # ---- SAM side: the test proteins (same order as in the bake-off) and the SAM train / validation proteins
    df, mats = b.load_data(rel)
    entries = df["Entry"].values
    pos = {e: i for i, e in enumerate(entries)}
    te = np.array([pos[e] for e in sam_sp["test"] if e in pos])
    X_te, y_te = mats["prott5"][te], df[LABELS].values.astype(np.float32)[te]
    if smoke:
        sel = np.random.RandomState(0).choice(len(te), min(150, len(te)), replace=False)
        X_te, y_te, te = X_te[sel], y_te[sel], te[sel]
        epochs, patience = 20, 5

    # ---- UP side: labels and ProtT5 vectors of this split's UP train / validation clusters
    lab = pd.read_csv(rel / "labels_UP.csv", usecols=["Entry", "in_benchmark", *LABELS], low_memory=False)
    lab = lab[lab.in_benchmark == True].drop_duplicates("Entry").set_index("Entry")                # noqa: E712
    tr_e, va_e = up_sp["train"], up_sp["validation"]
    if smoke:
        rng = np.random.RandomState(1)
        tr_e, va_e = list(rng.choice(tr_e, 3000, replace=False)), list(rng.choice(va_e, 600, replace=False))
    vecs = load_vectors(rel / "embeddings" / "embeddings_UP_ProtT5.h5", set(tr_e) | set(va_e))
    tr_ok, X_up_tr = stack(vecs, tr_e)
    va_ok, X_up_va = stack(vecs, va_e)
    missing = (len(tr_e) - len(tr_ok), len(va_e) - len(va_ok))
    if missing[0] or missing[1]:
        print(f"WARNING: {missing[0]} UP train and {missing[1]} UP validation proteins have no ProtT5 vector and are left out (embeddings/run_embeddings.sh, EMBED_WHAT=up-benchmark)", flush=True)
    y_up_tr, y_up_va = lab.loc[tr_ok, LABELS].values.astype(np.float32), lab.loc[va_ok, LABELS].values.astype(np.float32)
    print(f"iteration {it}: UP train {len(tr_ok)}, UP validation {len(va_ok)} (positives M/P/R/D train {y_up_tr.sum(0).astype(int).tolist()}, validation {y_up_va.sum(0).astype(int).tolist()}), "
          f"SAM test {len(te)} (positives {y_te.sum(0).astype(int).tolist()}) | epochs {epochs}, patience {patience}, threshold {b.THRESHOLD} | device {b.DEVICE} | smoke={smoke}", flush=True)

    overall, per_label, scores, curves, training = {}, {}, {"y_test": y_te, "test_entries": entries[te]}, [], []

    def train_arm(name, Xtr, ytr, Xva, yva, seed):
        t0 = time.time()
        model, hist, stop = b.train_early_stop(Xtr, ytr, Xva, yva, epochs, patience, b.LEARNING_RATE, seed)
        p_te = b.predict(model, X_te)
        metrics_rows(name, y_te, p_te, it, len(Xtr), len(Xva), len(X_te), overall, per_label)
        scores[f"{name}__test"] = p_te.astype(np.float32)
        ep = np.arange(1, epochs + 1)
        curves.append(pd.DataFrame({"iteration": it, "arm": name, "epoch": np.repeat(ep, len(LABELS)), "label": np.tile(LABELS, epochs),
                                    "train_loss": hist["train"].reshape(-1), "val_loss": hist["val"].reshape(-1)}))
        training.append({"iteration": it, "arm": name, "best_epoch": stop.best_epoch, "stop_epoch": stop.stop_epoch if stop.stop_epoch is not None else np.nan,
                         "best_val_total": stop.best, "val_total_at_end": float(hist["val"][-1].sum()), "epochs_trained": epochs, "n_train": len(Xtr), "n_val": len(Xva)})
        print(f"  {name:16s} MCC {overall[name]['mcc']:.3f} | train {len(Xtr)}, kept epoch {stop.best_epoch}, early stopping fired at {stop.stop_epoch} | {time.time() - t0:.0f}s", flush=True)

    if "SAM" in ARMS:
        z = np.load(bake, allow_pickle=True)
        if not smoke:
            assert np.array_equal(z["test_entries"], entries[te]), "the bake-off's test proteins differ from this split's SAM test set"
            sam_n = len(sam_sp["train"])
            metrics_rows("ProtT5_SAM", y_te, z["MLP_ProtT5__test"], it, sam_n, len(sam_sp["validation"]), len(te), overall, per_label)
            print(f"  {'ProtT5_SAM':16s} MCC {overall['ProtT5_SAM']['mcc']:.3f} | the bake-off's MLP_ProtT5 (SAM train {sam_n})", flush=True)
    if "UP" in ARMS:
        train_arm("ProtT5_UP", X_up_tr, y_up_tr, X_up_va, y_up_va, b.SEED + 100 * (it - 1) + 91)
    if "SAMclean" in ARMS:
        up = pd.read_csv(rel / "clusters" / "clusters_up_benchmark_id30.csv", usecols=["Entry", "cluster"]).set_index("Entry").cluster
        bad = set(up_sp["removed_clusters"])
        keep = lambda es: [e for e in es if e in pos and up.get(e, -1) not in bad]
        s_tr, s_va = keep(sam_sp["train"]), keep(sam_sp["validation"])
        n_drop = (len(sam_sp["train"]) - len(s_tr), len(sam_sp["validation"]) - len(s_va))
        if smoke:
            s_tr, s_va = s_tr[:600], s_va[:150]
        y_all = df[LABELS].values.astype(np.float32)
        itr, iva = np.array([pos[e] for e in s_tr]), np.array([pos[e] for e in s_va])
        print(f"  SAMclean: {n_drop[0]} SAM train and {n_drop[1]} validation proteins dropped (joint cluster with a test protein)", flush=True)
        train_arm("ProtT5_SAMclean", mats["prott5"][itr], y_all[itr], mats["prott5"][iva], y_all[iva], b.SEED + 100 * (it - 1) + 92)

    meta = {"split_mode": "cluster", "dataset": "full", "iteration": it}
    pd.DataFrame([{**meta, "classifier": c, **m} for c, m in overall.items()]).to_csv(f_det, index=False)
    pd.DataFrame([{**meta, "classifier": c, "label": lab_, **lm} for c, d in per_label.items() for lab_, lm in d.items()]).to_csv(f_pl, index=False)
    if curves:
        pd.concat(curves).to_csv(out_dir / f"iter{it:02d}_curves.csv.gz", index=False)
        pd.DataFrame(training).to_csv(out_dir / f"iter{it:02d}_training.csv", index=False)
    np.savez_compressed(out_dir / f"iter{it:02d}_scores.npz", **scores)
    print(f"wrote {f_det}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--iteration", type=int, required=True)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    import torch
    torch.set_num_threads(b.THREADS)
    run_iteration(a.iteration, a.smoke, a.force)


if __name__ == "__main__":
    main()
