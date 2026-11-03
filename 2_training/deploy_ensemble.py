#!/usr/bin/env python3
"""
Stage 2 -- the deployment model of the refresh: a ProtT5 MLP ensemble built from a K-fold CLUSTER partition (K = 6 members).

  python deploy_ensemble.py [--smoke]          # writes <release>/deployment/ ;  sbatch run_deploy_ensemble.sh

The design (chosen on the 2026_03 release, `ensemble_designs`, 16 held-out test sets): six members, EACH trained on a different 5/6 of the clusters and early-stopped on the other 1/6:
  * 6 seeds on the same data gave nothing over one model (MCC 0.9409 vs 0.9405); six members with different cluster splits gave 0.9558 (+0.015, better in 15 of 16 test sets, all four labels);
  * a K-fold PARTITION of the clusters was exactly as accurate as six independent random 15 % validation draws (0.9558 vs 0.9558), but it validates every protein exactly once, whereas six random
    draws leave 38 % of the proteins never validated -> here every protein gets one OUT-OF-FOLD prediction, from the one member that never saw its cluster;
  * the gain is saturated at about four to six members.
Data: labels_SAM.csv proteins with in_deployment (<= 2,000 aa), their clusters (cluster_deployment, mmseqs_clusters.py sam --subset deployment) and the ProtT5 vectors of
<release>/embeddings/embeddings_SAM_ProtT5.h5. The clusters are packed into K folds balanced by protein count AND by the positives of every label (greedy, best of 200 seeded orderings).
Every member is the bake-off's MLP (1024 -> 512 -> 256 -> 4, plain BCE, Adam 0.001, 2,000 epochs trained, best validation epoch kept, patience 80); settings: ../refresh_config.json.
Prediction = the mean of the six members' probabilities, THRESHOLD 0.5 per label (the benchmark's protocol); the member spread (SD) is the disagreement signal.

Writes <release>/deployment/ (the files 3_predictions/predict.py reads, so it is a drop-in: DEPLOY_DIR=<release>/deployment):
  ensemble_member_{0..K-1}.pt, model_config.json, thresholds.json (0.5 + the out-of-fold metrics), ad_reference.h5 (the training vectors) with the applicability-domain cutoff in the config
  (the 1st percentile of every training protein's nearest-neighbour cosine similarity to any other training protein, as in the archived model),
  oof_predictions.csv (Entry, cluster, fold, labels, oof_* = the honest out-of-fold probabilities, valfold_* = the deployed members' optimistic validation-fold probabilities), calibration.csv (reliability of oof_*), kept_epochs.csv,
  curves.csv.gz and deployment_train_val_loss.png (train / validation loss of every member over all epochs), report.txt.
HONEST OUT-OF-FOLD PREDICTIONS need a second set of six AUDIT members: a deployed member chose its stopping epoch on its own validation fold, so its predictions there (valfold_*) are optimistic (on
2026_03: MCC 0.97 against 0.94-0.95 on really unseen clusters, whatever the number of folds). Audit member j is trained like the deployed ones, but fold j is its TEST fold (never used for training or for
choosing the epoch) and fold (j + 1) mod K its validation fold; so every protein gets exactly one prediction (oof_*) from a model that neither trained on nor was selected on its cluster. Audit members train on
(K - 2) / K of the data, a little less than the deployed ones, so oof_* is conservative for the deployed ensemble. They are not part of the deployed model; oof_*, calibration.csv and the metrics in thresholds.json come from them.
"""
import argparse
import json
import os
import sys
from datetime import date
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch

import bakeoff as b

LABELS = b.LABELS
DEP = b.CFG.get("deployment", {})
N_MEMBERS = int(DEP.get("n_members", 6))
AD_PERCENTILE = float(DEP.get("ad_percentile", 1))


def pack_folds(clusters, y, k, seed, tries=200):
    """Assign whole clusters to k folds, balanced by protein count and by the positives of every label.
    -> (fold of every protein, [k, 1 + labels] load of every fold, the worst relative deviation from the ideal load)."""
    ids, inv = np.unique(clusters, return_inverse=True)
    size = np.bincount(inv).astype(float)
    lab = np.zeros((len(ids), y.shape[1]))
    np.add.at(lab, inv, y)
    metric = np.column_stack([size, lab])
    target = np.maximum(metric.sum(0) / k, 1.0)
    rng = np.random.default_rng(np.random.SeedSequence([seed, 5]))
    best = None
    for _ in range(tries):
        order = np.lexsort((rng.random(len(ids)), -size))                # big clusters first, random among equal sizes
        load, assign = np.zeros((k, metric.shape[1])), np.empty(len(ids), int)
        for c in order:
            f = int(np.argmin((((load + metric[c]) / target) ** 2).sum(1)))      # the fold whose normalised load grows least
            assign[c] = f
            load[f] += metric[c]
        worst = float(np.abs(load / target - 1).max())
        if best is None or worst < best[0]:
            best = (worst, assign.copy(), load.copy())
    return best[1][inv], best[2], best[0]


def calibration_table(y, p, bins=10):
    rows = []
    edges = np.linspace(0, 1, bins + 1)
    for j, lab in enumerate(LABELS):
        which = np.clip(np.digitize(p[:, j], edges[1:-1]), 0, bins - 1)
        for k in range(bins):
            m = which == k
            if m.sum():
                rows.append({"label": lab, "bin": f"{edges[k]:.1f}-{edges[k + 1]:.1f}", "n": int(m.sum()), "mean_predicted": float(p[m, j].mean()), "observed_positive_rate": float(y[m, j].mean())})
    return pd.DataFrame(rows)


def ad_cutoff_from(X, percentile):
    """1st percentile (default) of the nearest-neighbour cosine similarity of every training vector to any OTHER training vector."""
    Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-8)
    t = torch.tensor(Xn, device=b.DEVICE)
    nn_sim = np.empty(len(Xn), dtype=np.float32)
    for i in range(0, len(Xn), 2048):
        s = t[i:i + 2048] @ t.T
        s[torch.arange(s.shape[0]), torch.arange(i, i + s.shape[0])] = -2.0                  # not itself
        nn_sim[i:i + 2048] = s.max(1).values.cpu().numpy()
    return float(np.percentile(nn_sim, percentile)), nn_sim


def load_ensemble(model_dir, device=None):
    """The members of a deployment folder -> (config, [MultilabelNN in eval mode])."""
    model_dir = Path(model_dir)
    cfg = json.loads((model_dir / "model_config.json").read_text())
    models = []
    for j in range(cfg["n_ensemble"]):
        m = b.MultilabelNN(cfg["input_dim"], len(cfg["labels"])).to(device or b.DEVICE)
        m.load_state_dict(torch.load(model_dir / f"ensemble_member_{j}.pt", map_location=device or b.DEVICE))
        models.append(m.eval())
    return cfg, models


def ensemble_predict(models, X):
    """-> (mean probability [n, labels], member SD [n, labels])."""
    P = np.stack([b.predict(m, X) for m in models])
    return P.mean(0), P.std(0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smoke", action="store_true", help="20 epochs, 1,500 proteins; writes to <release>/deployment/_smoke")
    a = ap.parse_args()
    rel = b.release_dir()
    out = Path(os.environ.get("DEPLOY_OUT", rel / "deployment")) / ("_smoke" if a.smoke else "")
    out.mkdir(parents=True, exist_ok=True)
    lab = pd.read_csv(rel / "labels_SAM.csv", usecols=["Entry", "in_deployment", "cluster_deployment", *LABELS], low_memory=False)
    lab = lab[lab.in_deployment == True].drop_duplicates("Entry").reset_index(drop=True)             # noqa: E712
    if lab.cluster_deployment.isna().any():
        sys.exit("some deployment proteins have no cluster: run 1_clustering/run_mmseqs_clusters.sh (sam-deployment) first")
    vec_file = rel / "embeddings" / "embeddings_SAM_ProtT5.h5"
    vecs = b.load_source(vec_file, "flat")
    n0 = len(lab)
    lab = lab[lab.Entry.isin(vecs)].reset_index(drop=True)
    if len(lab) < n0:
        print(f"WARNING: {n0 - len(lab)} of {n0} deployment proteins have no ProtT5 vector and are left out (embeddings/run_embeddings.sh)", flush=True)
    epochs, patience = b.EPOCHS, b.PATIENCE
    if a.smoke:
        lab = lab.sample(1500, random_state=0).reset_index(drop=True)
        epochs, patience = 20, 5
    X = np.vstack([vecs[e] for e in lab.Entry]).astype(np.float32)
    y = lab[LABELS].values.astype(np.float32)
    fold, load, worst = pack_folds(lab.cluster_deployment.values, y, N_MEMBERS, b.SEED)
    assert len(np.unique(fold)) == N_MEMBERS
    lines = []
    say = lambda s="": (print(s, flush=True), lines.append(s))
    say(f"deployment set: {len(lab)} proteins in {lab.cluster_deployment.nunique()} clusters, ProtT5 {X.shape[1]}-D | {N_MEMBERS} members, one per fold | epochs {epochs}, patience {patience} | device {b.DEVICE}")
    say(f"folds (balanced by size and labels; worst deviation from the ideal load {worst:.1%}):  proteins (Metabolite / Protein / RNA / DNA positives) of the validation fold of each member")
    for j in range(N_MEMBERS):
        m = fold == j
        say(f"   member {j}: {int(m.sum()):5d} ({int(y[m, 0].sum())} / {int(y[m, 1].sum())} / {int(y[m, 2].sum())} / {int(y[m, 3].sum())})")

    members, oof, kept, curves = [], np.zeros_like(y), [], []
    for j in range(N_MEMBERS):
        va, tr = fold == j, fold != j
        model, hist, stop = b.train_early_stop(X[tr], y[tr], X[va], y[va], epochs, patience, b.LEARNING_RATE, b.SEED + 1000 + j)
        oof[va] = b.predict(model, X[va])
        members.append(model)
        torch.save({k: v.detach().cpu() for k, v in model.state_dict().items()}, out / f"ensemble_member_{j}.pt")
        kept.append({"member": j, "best_epoch": stop.best_epoch, "stop_epoch": stop.stop_epoch if stop.stop_epoch is not None else np.nan, "best_val_total": stop.best,
                     "val_total_at_end": float(hist["val"][-1].sum()), "n_train": int(tr.sum()), "n_val": int(va.sum())})
        ep = np.arange(1, epochs + 1)
        curves.append(pd.DataFrame({"member": j, "epoch": np.repeat(ep, len(LABELS)), "label": np.tile(LABELS, epochs), "train_loss": hist["train"].reshape(-1), "val_loss": hist["val"].reshape(-1)}))
        say(f"member {j}: trained on {int(tr.sum())}, validated on {int(va.sum())}; kept epoch {stop.best_epoch}, early stopping fired at {stop.stop_epoch}")

    vm = b.all_metrics(y, oof, 0.5)
    say(f"\nvalidation-fold predictions of the deployed members (OPTIMISTIC: the stopping epoch was chosen on these very proteins): MCC {vm['mcc']:.4f}")
    valfold = oof.copy()
    oof = np.zeros_like(y)
    audit_kept = []
    for j in range(N_MEMBERS):                                                     # audit members: test fold j, validation fold j+1, train the rest
        te, va = fold == j, fold == (j + 1) % N_MEMBERS
        tr = ~(te | va)
        model, _, stop = b.train_early_stop(X[tr], y[tr], X[va], y[va], epochs, patience, b.LEARNING_RATE, b.SEED + 2000 + j)
        oof[te] = b.predict(model, X[te])
        audit_kept.append(stop.best_epoch)
    om, opl = b.all_metrics(y, oof, 0.5), b.per_label_metrics(y, oof, 0.5)
    per_fold = [b.all_metrics(y[fold == j], oof[fold == j], 0.5)["mcc"] for j in range(N_MEMBERS)]
    say(f"HONEST out-of-fold (audit members: fold j never used for training or for the epoch; fixed 0.5): pooled MCC {om['mcc']:.4f}, mean over folds {np.mean(per_fold):.4f} (SD {np.std(per_fold, ddof=1):.3f}), "
        f"AUPRC {om['auprc']:.4f}, AUROC {om['auroc']:.4f}, F1 {om['f1']:.4f}; the optimism of the validation-fold numbers: {vm['mcc'] - om['mcc']:+.4f}")
    say("   per label MCC: " + ", ".join(f"{l} {opl[l]['mcc']:.3f}" for l in LABELS) + f" | audit members' kept epochs {audit_kept}")
    ad_cut, nn_sim = ad_cutoff_from(X, AD_PERCENTILE)
    say(f"applicability domain: nearest-neighbour cosine among the training proteins, {AD_PERCENTILE:g}th percentile = {ad_cut:.4f} (median {np.median(nn_sim):.4f})")

    with h5py.File(out / "ad_reference.h5", "w") as f:
        for e, v in zip(lab.Entry, X):
            f.create_dataset(e, data=v)
    (out / "thresholds.json").write_text(json.dumps({"thresholds": {l: 0.5 for l in LABELS}, "rule": "fixed 0.5 on the mean probability of the members (the benchmark's protocol)",
                                                     "oof_metrics": {"note": "honest: audit members whose test fold was never used for training or epoch selection", "mcc": om["mcc"], "mcc_mean_over_folds": float(np.mean(per_fold)), "auprc": om["auprc"], "auroc": om["auroc"],
                                                                  "per_label_mcc": {l: opl[l]["mcc"] for l in LABELS}, "optimistic_validation_fold_mcc": vm["mcc"]}}, indent=1))
    (out / "model_config.json").write_text(json.dumps({
        "embedding": "ProtT5", "embedding_h5": vec_file.name, "input_dim": int(X.shape[1]), "hidden": [512, 256], "labels": LABELS, "loss": "bce (all labels)", "n_ensemble": N_MEMBERS,
        "design": f"{N_MEMBERS}-fold cluster partition: member j validated on fold j, trained on the other folds; best validation epoch kept (patience {patience}, {epochs} epochs trained); honest out-of-fold predictions from separate audit members",
        "kept_epochs": [k["best_epoch"] for k in kept], "epochs": epochs, "patience": patience, "learning_rate": b.LEARNING_RATE, "trained_on": int(len(lab)),
        "ad_reference_h5": "ad_reference.h5", "ad_reference_path": "ad_reference.h5", "ad_cutoff": ad_cut, "ad_percentile": AD_PERCENTILE, "built": str(date.today()), "release": rel.name,
        "refresh_config": b.CFG}, indent=1))
    pd.DataFrame({"Entry": lab.Entry, "cluster": lab.cluster_deployment, "fold": fold, **{l: lab[l] for l in LABELS}, **{f"oof_{l}": oof[:, j].round(5) for j, l in enumerate(LABELS)}, **{f"valfold_{l}": valfold[:, j].round(5) for j, l in enumerate(LABELS)}}).to_csv(out / "oof_predictions.csv", index=False)
    calibration_table(y, oof).to_csv(out / "calibration.csv", index=False)
    pd.DataFrame(kept).to_csv(out / "kept_epochs.csv", index=False)
    cur = pd.concat(curves)
    cur.to_csv(out / "curves.csv.gz", index=False)
    (out / "report.txt").write_text("\n".join(lines) + "\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 3, figsize=(15, 7), squeeze=False)
    for ax, j in zip(axs.ravel(), range(N_MEMBERS)):
        c = cur[cur.member == j].groupby("epoch")[["train_loss", "val_loss"]].sum()
        ax.plot(c.index, c.train_loss, color="#1f77b4", label="train"); ax.plot(c.index, c.val_loss, color="#d62728", label="validation")
        ax.axvline(kept[j]["best_epoch"], color="k", ls="--", lw=1); ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("epoch"); ax.set_ylabel("total loss")
        ax.set_title(f"member {j}: kept epoch {kept[j]['best_epoch']}", fontsize=9)
    axs[0, 0].legend(fontsize=8)
    fig.suptitle("Deployment ensemble: train / validation loss of the six members (dashed = the epoch kept)", fontsize=10)
    fig.tight_layout(); fig.savefig(out / "deployment_train_val_loss.png", dpi=120); plt.close(fig)
    say(f"\nwrote {out}/ ({N_MEMBERS} members, thresholds.json, model_config.json, ad_reference.h5, oof_predictions.csv, calibration.csv, kept_epochs.csv, curves.csv.gz, figure, report.txt)")


if __name__ == "__main__":
    main()
