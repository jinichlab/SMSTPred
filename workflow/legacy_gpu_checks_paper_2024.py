"""Checks for the GPU re-runs (tier B steps of `make verify-gpu`).

A retrained neural network does not reproduce to the last digit across PyTorch / CUDA / GPU generations, so these steps compare outcomes with
stated tolerances instead of comparing files. The tolerances are set from measured noise: the committed models were trained in an earlier
environment, and re-training them with environment.yml on an RTX 6000 Ada moved thresholds by up to 0.008, out-of-fold scores by up to 0.012,
bake-off MCCs by up to 0.04 per iteration (0.005 per method) and 0.2% of the class calls. Training is deterministic within one environment
(two runs gave identical results).

Each check takes (work, repo): the scratch copy that was run and the repository holding the committed results, and returns a list of
(what, ok, detail).
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

LABELS = ["Metabolite", "Protein", "RNA", "DNA"]


def bakeoff(work, repo):
    key = ["split_mode", "dataset", "iteration", "classifier"]
    ref = pd.read_csv(repo / "2_training/sam_classifier_comparison_output/sam_classifier_comparison_detailed.csv")
    new = pd.read_csv(work / "2_training/sam_classifier_comparison_output_mlp_only/sam_classifier_comparison_detailed.csv")
    m = ref.merge(new, on=key, suffixes=("_ref", "_new"))
    res = [("every MLP row of the re-run has a committed counterpart", len(m) == len(new) > 0, f"{len(m)} of {len(new)} rows")]
    d = (m.mcc_new - m.mcc_ref).abs()
    res.append(("max absolute MCC difference per iteration at most 0.06", bool(d.max() <= 0.06), f"max {d.max():.3f}, mean {d.mean():.4f}"))
    g = m.groupby("classifier")[["mcc_ref", "mcc_new"]].mean()
    dg = (g.mcc_new - g.mcc_ref).abs()
    res.append(("mean MCC per method within 0.01 of the committed run", bool((dg <= 0.01).all()),
                "; ".join(f"{c} {g.mcc_new[c]:.4f} vs {g.mcc_ref[c]:.4f}" for c in g.index)))
    return res


def deployment(work, repo):
    a, b = repo / "2_training/deployment_model", work / "2_training/deployment_model"
    ta, tb = json.loads((a / "thresholds.json").read_text()), json.loads((b / "thresholds.json").read_text())
    oa, ob = pd.read_csv(a / "oof_predictions.csv"), pd.read_csv(b / "oof_predictions.csv")
    res = [("same proteins, clusters and folds", bool((oa[["Entry", "cluster", "fold"]] == ob[["Entry", "cluster", "fold"]]).all().all()), f"{len(oa)} proteins")]
    dt = {k: abs(ta["thresholds"][k] - tb["thresholds"][k]) for k in LABELS}
    res.append(("per-label thresholds within 0.02", max(dt.values()) <= 0.02, ", ".join(f"{k} {tb['thresholds'][k]:.3f} vs {ta['thresholds'][k]:.3f}" for k in LABELS)))
    dm = {k: abs(ta["oof_metrics"][k]["mcc"] - tb["oof_metrics"][k]["mcc"]) for k in LABELS}
    res.append(("per-label out-of-fold MCC within 0.01", max(dm.values()) <= 0.01, ", ".join(f"{k} {tb['oof_metrics'][k]['mcc']:.3f} vs {ta['oof_metrics'][k]['mcc']:.3f}" for k in LABELS)))
    ds = max(float(np.abs(oa[f"score_{k}"] - ob[f"score_{k}"]).max()) for k in LABELS)
    res.append(("out-of-fold scores within 0.03", ds <= 0.03, f"max difference {ds:.4f}"))
    return res


def predictions(work, repo):
    a, b = repo / "3_predictions/predictions_output", work / "3_predictions/predictions_output"
    fa, fb = sorted(p.name for p in a.glob("*.csv")), sorted(p.name for p in b.glob("*.csv"))
    res = [("same organism files", fa == fb, f"{len(fb)} files")]
    ra, rb = [], []
    for f in fa:
        x, y = pd.read_csv(a / f), pd.read_csv(b / f)
        assert list(x.Entry) == list(y.Entry), f
        ra.append(x), rb.append(y)
    x, y = pd.concat(ra, ignore_index=True), pd.concat(rb, ignore_index=True)
    calls = [f"call_{k}" for k in LABELS]
    same = float((x[calls].values == y[calls].values).all(1).mean())
    res.append(("identical class calls for at least 99.5% of proteins", same >= 0.995, f"{100 * same:.2f}% of {len(x):,} proteins"))
    ds = max(float(np.abs(x[f"score_{k}"] - y[f"score_{k}"]).max()) for k in LABELS)
    res.append(("scores within 0.05", ds <= 0.05, f"max difference {ds:.4f}"))
    res.append(("applicability-domain flags identical", bool((x.in_domain == y.in_domain).all()), ""))
    return res


def loss_curves(work, repo):
    a, b = repo / "2_training/train_val_loss_output", work / "2_training/train_val_loss_output"
    files = sorted(p.name for p in a.glob("*__curves.csv.gz") if not p.name.startswith("sam_balanced"))   # the balanced sweeps are not in the manuscript
    got = sorted(p.name for p in b.glob("*__curves.csv.gz"))
    miss = [f for f in files if f not in got]
    res = [("a curves file for every committed sweep", not miss, f"{len(files) - len(miss)} of {len(files)} sweeps" + (f"; missing {miss}" if miss else ""))]
    k = ["dataset", "embedding", "loss_type", "iteration", "epoch", "label"]
    dv, dm, fin = [], [], []
    for f in got:
        if f not in files:
            continue
        x, y = pd.read_csv(a / f), pd.read_csv(b / f)
        m = x.merge(y, on=k, suffixes=("_ref", "_new"))
        dv.append((m.val_loss_new - m.val_loss_ref).abs())
        last = m[m.epoch == m.epoch.max()]
        fin.append((last.val_loss_new - last.val_loss_ref).abs())
        c = f.replace("__curves", "__ckpt_mcc")
        if (b / c).exists():
            p, q = pd.read_csv(a / c), pd.read_csv(b / c)
            mm = p.merge(q, on=[z for z in p.columns if z != "mcc"], suffixes=("_ref", "_new"))
            dm.append((mm.mcc_new - mm.mcc_ref).abs())
    dv, fin, dm = pd.concat(dv), pd.concat(fin), pd.concat(dm)
    res.append(("median absolute validation-loss difference over all epochs at most 0.01", float(dv.median()) <= 0.01, f"median {dv.median():.4f}, 95th percentile {dv.quantile(.95):.4f}"))
    res.append(("median absolute final-epoch validation-loss difference at most 0.01", float(fin.median()) <= 0.01, f"median {fin.median():.4f}, max {fin.max():.4f}"))
    res.append(("median absolute checkpoint-MCC difference at most 0.02", float(dm.median()) <= 0.02, f"median {dm.median():.4f}, 95th percentile {dm.quantile(.95):.4f}"))
    return res


CHECKS = {"bakeoff": bakeoff, "deployment": deployment, "predictions": predictions, "loss_curves": loss_curves}
