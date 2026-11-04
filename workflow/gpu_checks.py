"""Checks for the GPU re-runs (tier B steps of `make verify-gpu`).

A retrained neural network does not reproduce to the last digit across PyTorch / CUDA / GPU generations, so these steps compare outcomes with stated tolerances
instead of comparing files. Each check takes (work_release, ref_release): the scratch copy of the release that was run and the archived release, and returns a list
of (what, ok, detail). The tolerances below are set from the run-to-run and environment-to-environment spread measured on the 2024-era model (see
workflow/legacy_gpu_checks_paper_2024.py: thresholds 0.008, out-of-fold scores 0.012, per-label MCC about 0.01) and widened for six members; they are provisional until
several GPU generations have been compared, and a run that exceeds them is reported with the measured differences rather than silently relaxed.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

LABELS = ["Metabolite", "Protein", "RNA", "DNA"]


def deployment(work, ref):
    a, b = Path(ref) / "deployment", Path(work) / "deployment"
    ca, cb = json.loads((a / "model_config.json").read_text()), json.loads((b / "model_config.json").read_text())
    ta, tb = json.loads((a / "thresholds.json").read_text()), json.loads((b / "thresholds.json").read_text())
    oa, ob = pd.read_csv(a / "oof_predictions.csv"), pd.read_csv(b / "oof_predictions.csv")
    res = [("same proteins, clusters and folds", bool((oa[["Entry", "cluster", "fold"]] == ob[["Entry", "cluster", "fold"]]).all().all()), f"{len(oa)} proteins")]
    res.append(("six ensemble members written", all((b / f"ensemble_member_{i}.pt").exists() for i in range(6)), ""))
    res.append(("every label called at the same fixed threshold as archived", ta["thresholds"] == tb["thresholds"], str(tb["thresholds"])))
    ma, mb = ta["oof_metrics"], tb["oof_metrics"]
    res.append(("overall out-of-fold MCC within 0.01", abs(ma["mcc"] - mb["mcc"]) <= 0.01, f"{mb['mcc']:.4f} vs {ma['mcc']:.4f}"))
    dm = {k: abs(ma["per_label_mcc"][k] - mb["per_label_mcc"][k]) for k in LABELS}
    res.append(("per-label out-of-fold MCC within 0.02", max(dm.values()) <= 0.02, ", ".join(f"{k} {mb['per_label_mcc'][k]:.3f} vs {ma['per_label_mcc'][k]:.3f}" for k in LABELS)))
    ds = max(float(np.abs(oa[f"oof_{k}"] - ob[f"oof_{k}"]).max()) for k in LABELS if f"oof_{k}" in oa)
    res.append(("out-of-fold scores within 0.1 (max over proteins)", ds <= 0.1, f"max difference {ds:.4f}"))
    dc = float(np.mean([np.abs(oa[f"oof_{k}"] - ob[f"oof_{k}"]).mean() for k in LABELS if f"oof_{k}" in oa]))
    res.append(("mean out-of-fold score difference within 0.01", dc <= 0.01, f"mean difference {dc:.4f}"))
    if "ad_cutoff" in ca and "ad_cutoff" in cb:
        res.append(("applicability-domain cutoff within 0.01", abs(ca["ad_cutoff"] - cb["ad_cutoff"]) <= 0.01, f"{cb['ad_cutoff']:.4f} vs {ca['ad_cutoff']:.4f}"))
    return res


CHECKS = {"deployment": deployment}
