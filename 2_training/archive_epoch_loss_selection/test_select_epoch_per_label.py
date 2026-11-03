#!/usr/bin/env python3
"""Tests of select_epoch_per_label.py on synthetic curves: every (embedding, label) gets ITS OWN minimum, a late dip after overfitting is ignored, one split is enough."""
from types import SimpleNamespace

import numpy as np
import pandas as pd

import select_epoch_per_label as s

A = SimpleNamespace(tol=0.03, patience=100, step=10, smooth=7, max_search=3000)


def curves(minima, n_split=3, late_dip=False, seed=0):
    """BCE validation loss per (embedding, label): falls to its minimum EXACTLY at the given epoch, then rises; a deep dip near epoch 2500 if late_dip."""
    rng = np.random.default_rng(seed)
    ep = np.arange(1, 3001)
    rows = []
    for (emb, lab), m in minima.items():
        for it in range(1, n_split + 1):
            v = 0.05 + 0.5 * (np.maximum(m - ep, 0) / m) ** 2 + 0.00005 * np.maximum(ep - m, 0) + rng.normal(0, 0.0004, len(ep))   # exactly at epoch m
            if late_dip:
                v[2480:2520] = 0.0
            rows.append(pd.DataFrame({"embedding": emb, "loss_type": "bce", "iteration": it, "epoch": ep, "label": lab, "val_loss": v, "train_loss": v * 0.5}))
    return pd.concat(rows, ignore_index=True)


def test_each_cell_gets_its_own_minimum():
    minima = {("E1", "Protein"): 150, ("E1", "DNA"): 400, ("E2", "Protein"): 250, ("E2", "DNA"): 90}
    t = s.select(s.label_cells(curves(minima), "bce"), A, 20)
    got = {(r.embedding, r.label): r.epoch for r in t.itertuples()}
    for k, m in minima.items():
        assert abs(got[k] - m) <= 60, (k, got[k], m)          # the smoothed minimum sits within a few tens of epochs of the true one
    assert got[("E1", "DNA")] > got[("E1", "Protein")] > 0 and got[("E2", "DNA")] < got[("E2", "Protein")]   # and the order of the cells is kept
    assert t.onset_found.all()


def test_late_dip_is_ignored():
    t = s.select(s.label_cells(curves({("E1", "Protein"): 200}, late_dip=True), "bce"), A, 5)
    assert abs(t.epoch.iloc[0] - 200) <= 60 and t.epoch.iloc[0] < 1000   # the dip to 0 at epoch 2500 is after overfitting set in: never looked at


def test_one_split_is_enough():
    t = s.select(s.label_cells(curves({("E1", "RNA"): 300}, n_split=1), "bce"), A, 5)
    assert abs(t.epoch.iloc[0] - 300) <= 80


def test_no_bootstrap_is_allowed():
    t = s.select(s.label_cells(curves({("E1", "RNA"): 300}), "bce"), A, 0)
    assert t.ci95_lo.iloc[0] == t.ci95_hi.iloc[0] == t.epoch.iloc[0]                  # --bootstrap 0: the interval collapses to the epoch instead of crashing


if __name__ == "__main__":
    for t in (test_each_cell_gets_its_own_minimum, test_late_dip_is_ignored, test_one_split_is_enough, test_no_bootstrap_is_allowed):
        t()
        print("ok", t.__name__)
