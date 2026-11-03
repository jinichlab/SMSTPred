#!/usr/bin/env python3
"""The fast threshold search of train_val_loss_curves.py must return exactly the thresholds of the original loop (every candidate scored with sklearn's matthews_corrcoef,
first maximum kept). Hard cases: ties from repeated probabilities, a very rare class, a single class, perfect and useless predictors. Run: python test_threshold_search.py"""
import time

import numpy as np
from sklearn.metrics import matthews_corrcoef

import train_val_loss_curves as t


def original(y_true, probs):
    thr = np.full(len(t.LABELS), 0.5)
    for i in range(len(t.LABELS)):
        if len(np.unique(y_true[:, i])) < 2:
            continue
        best_mcc, best_t = -2.0, 0.5
        for th in np.unique(np.percentile(probs[:, i], np.linspace(0, 100, 200))):
            m = matthews_corrcoef(y_true[:, i], (probs[:, i] >= th).astype(int))
            if m > best_mcc:
                best_mcc, best_t = m, th
        thr[i] = best_t
    return thr


def main():
    rng = np.random.default_rng(1)
    n_bad = 0
    n_cases = 300
    for k in range(n_cases):
        n = int(rng.integers(30, 3000))
        prev = rng.choice([0.0, 0.003, 0.015, 0.2, 0.5, 0.9, 1.0], size=4)
        y = (rng.random((n, 4)) < prev).astype(np.float32)
        noise = rng.choice([0.0, 0.05, 0.4, 1.0])
        p = np.clip(y * rng.random() + rng.random((n, 4)) * noise, 0, 1).astype(np.float32)
        if k % 3 == 0:
            p = np.round(p, 1)                 # many tied probabilities: the tie-break matters
        if k % 7 == 0:
            p[:, 0] = 0.5                      # a constant predictor
        a, b = t.find_mcc_optimal_thresholds(y, p), original(y, p)
        if not np.array_equal(a, b):
            n_bad += 1
            print("DIFFERENT", k, a, b)
    print(f"{n_cases} cases, {n_bad} differ")
    assert n_bad == 0
    y = (rng.random((6800, 4)) < [0.2, 0.15, 0.65, 0.015]).astype(np.float32)
    p = np.clip(y * 0.7 + rng.random((6800, 4)) * 0.4, 0, 1).astype(np.float32)
    t0 = time.time(); [t.find_mcc_optimal_thresholds(y, p) for _ in range(5)]; fast = (time.time() - t0) / 5
    t0 = time.time(); [original(y, p) for _ in range(2)]; slow = (time.time() - t0) / 2
    print(f"6,800 proteins: {fast * 1000:.0f} ms vs {slow * 1000:.0f} ms ({slow / fast:.0f}x faster)")


if __name__ == "__main__":
    main()
