#!/usr/bin/env python3
"""
Stage 2 -- MLP vs classical machine learning on the SAME ProtT5 embeddings (Fig S4 of the refreshed release): ONE split per run, the 16 splits of the bake-off.

  python classical_ml.py --iteration 3            # appends iter03_classical_{detailed,per_label}.csv to <release>/benchmark/iterations/
  python classical_ml.py --iteration 1 --smoke    # tiny end-to-end test
  CLASSICAL_METHODS="XGBoost,Logistic Regression" python classical_ml.py --iteration 1
  sbatch run_classical.sh (16-task CPU array) ; then sbatch run_bakeoff_report.sh (it merges the classical rows with the bake-off and draws Fig S4)

Seven classifiers, the settings of the archived paper's Fig S4 (scikit-learn defaults except those listed): XGBoost, Random Forest (100 trees), Logistic Regression (max_iter 1000),
SVM linear / RBF / polynomial (degree 3) / sigmoid (probability=True: Platt scaling by internal 5-fold CV). One model per label (MultiOutputClassifier; SVMs: one per label in parallel).
Features: the raw ProtT5 embeddings (1024-D, no scaling). Protocol identical to the MLP's: train on the SAME train clusters, the SAME test clusters, a fixed 0.5 threshold on the
predicted probability. No hyperparameter is tuned and the validation clusters are not used (these models have no early-stopping epoch; they keep their defaults, as in the archived figure) --
the MLP's only tuned quantity is its stopping epoch. The MLP_ProtT5 row it is compared with is the bake-off's. Reads the release folder and ../refresh_config.json (bakeoff.py).
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

import bakeoff as b

LABELS = b.LABELS
CLASSIFIERS = ["XGBoost", "Random Forest", "Logistic Regression", "SVM_Linear", "SVM_RBF", "SVM_Polynomial", "SVM_Sigmoidal"]


def templates(n_jobs):
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.svm import SVC
    from xgboost import XGBClassifier
    return {"XGBoost": XGBClassifier(eval_metric="logloss", random_state=42, n_jobs=n_jobs),
            "Random Forest": RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=n_jobs),
            "Logistic Regression": LogisticRegression(max_iter=1000, random_state=42, n_jobs=n_jobs),
            "SVM_Linear": SVC(kernel="linear", probability=True, random_state=42),
            "SVM_RBF": SVC(kernel="rbf", probability=True, random_state=42),
            "SVM_Polynomial": SVC(kernel="poly", degree=3, probability=True, random_state=42),
            "SVM_Sigmoidal": SVC(kernel="sigmoid", probability=True, random_state=42)}


def fit_predict(name, template, X_train, y_train, X_test, n_jobs):
    """Probability of the positive class for every test protein and label; a label with one class in the train set gets 0."""
    from joblib import Parallel, delayed
    from sklearn.base import clone
    from sklearn.multioutput import MultiOutputClassifier
    n_labels = y_train.shape[1]
    probs = np.zeros((X_test.shape[0], n_labels))
    if name.startswith("SVM"):
        def one(j):
            if len(np.unique(y_train[:, j])) < 2:
                return None
            return clone(template).fit(X_train, y_train[:, j])
        for j, m in enumerate(Parallel(n_jobs=min(n_jobs, n_labels))(delayed(one)(j) for j in range(n_labels))):
            if m is not None:
                probs[:, j] = m.predict_proba(X_test)[:, 1]
    else:
        model = MultiOutputClassifier(clone(template)).fit(X_train, y_train)
        probs = np.column_stack([e.predict_proba(X_test)[:, 1] for e in model.estimators_])
    return probs


def run_iteration(it, smoke, force):
    rel = b.release_dir()
    out_dir = Path(os.environ.get("BAKEOFF_OUT", rel / "benchmark")) / ("_smoke" if smoke else "") / "iterations"
    out_dir.mkdir(parents=True, exist_ok=True)
    f_det, f_pl = out_dir / f"iter{it:02d}_classical_detailed.csv", out_dir / f"iter{it:02d}_classical_per_label.csv"
    if f_det.exists() and not force:
        print(f"iteration {it}: {f_det.name} exists, skipping (use --force to redo)")
        return
    methods = [m.strip() for m in os.environ.get("CLASSICAL_METHODS", ",".join(CLASSIFIERS)).split(",")]
    assert all(m in CLASSIFIERS for m in methods), f"unknown classifier in {methods}; choose from {CLASSIFIERS}"
    sp_cfg = b.CFG["splits"]
    sp_file = rel / "clusters" / f"splits_tvt_sam_benchmark_id30_t{round(sp_cfg['test_frac'] * 100)}_v{round(sp_cfg['val_frac'] * 100)}.json"
    sp = json.loads(sp_file.read_text())["iterations"][it - 1]
    df, mats = b.load_data(rel)
    pos = {e: i for i, e in enumerate(df["Entry"].values)}
    y = df[LABELS].values.astype(np.float32)
    tr, te = (np.array([pos[e] for e in sp[part] if e in pos]) for part in ("train", "test"))
    if smoke:
        rng = np.random.RandomState(0)
        tr, te = rng.choice(tr, 600, replace=False), rng.choice(te, 150, replace=False)
    X_tr, X_te, y_tr, y_te = mats["prott5"][tr], mats["prott5"][te], y[tr], y[te]
    n_jobs = b.THREADS
    print(f"iteration {it}: ProtT5 features, train {len(tr)}, test {len(te)} | positives M/P/R/D train {y_tr.sum(0).astype(int).tolist()} test {y_te.sum(0).astype(int).tolist()} | threshold {b.THRESHOLD} | jobs {n_jobs} | smoke={smoke}", flush=True)
    tmpl = templates(n_jobs)
    meta = {"split_mode": "cluster", "dataset": "full", "iteration": it}
    det, pl, scores = [], [], {"y_test": y_te, "test_entries": df["Entry"].values[te]}
    for name in methods:
        t0 = time.time()
        p = fit_predict(name, tmpl[name], X_tr, y_tr, X_te, n_jobs)
        det.append({**meta, "classifier": name, "n_train": len(tr), "n_val": 0, "n_test": len(te), **b.all_metrics(y_te, p, b.THRESHOLD)})
        pl += [{**meta, "classifier": name, "label": lab, **m} for lab, m in b.per_label_metrics(y_te, p, b.THRESHOLD).items()]
        scores[f"{name}__test"] = p.astype(np.float32)
        print(f"  {name:20s} MCC {det[-1]['mcc']:.3f} | {time.time() - t0:.0f}s", flush=True)
    pd.DataFrame(det).to_csv(f_det, index=False)
    pd.DataFrame(pl).to_csv(f_pl, index=False)
    np.savez_compressed(out_dir / f"iter{it:02d}_classical_scores.npz", **scores)
    print(f"wrote {f_det}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--iteration", type=int, required=True)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    run_iteration(a.iteration, a.smoke, a.force)


if __name__ == "__main__":
    main()
