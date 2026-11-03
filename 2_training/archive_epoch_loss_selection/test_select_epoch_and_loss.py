#!/usr/bin/env python3
"""Tests of select_epoch_and_loss.py on synthetic data: the overfitting guard (a late dip or spikes must not be picked), the checkpoint snapping, the loss scopes and the
selection optimism. Run: python test_select_epoch_and_loss.py   (or pytest)"""
import argparse

import numpy as np
import pandas as pd

import select_epoch_and_loss as s

A = argparse.Namespace(tol=0.03, patience=100, max_search=3000, smooth=7, step=10)
LOG = np.array(list(range(1, 1001)) + list(range(1010, 20001, 10)), dtype=float)


def u_curve(rng, late_dip=True, spikes=True, noise=0.002):
    """Validation loss that falls fast, bottoms out near epoch 183 and then rises (overfitting); optionally with a deep dip and spikes late in training."""
    v = 0.25 + 0.75 * np.exp(-LOG / 50) + 0.35 * (1 - np.exp(-LOG / 700)) + rng.normal(0, noise, len(LOG))
    if late_dip:
        v[(LOG >= 15000) & (LOG <= 15100)] = 0.2        # far below the true minimum (about 0.35)
    if spikes:
        v[rng.choice(np.flatnonzero(LOG > 5000), 60, replace=False)] *= 2.5
    return v


def test_guard_ignores_late_dips_and_spikes():
    rng = np.random.default_rng(0)
    mat = np.stack([u_curve(rng) for _ in range(6)], axis=1)
    cols = np.arange(6)
    ep, found = s.cell_epoch(LOG, mat, cols, A)
    assert found and 150 <= ep <= 230, ep                       # the true minimum, near 183
    grid, y = s.smooth_curve(LOG, mat.mean(axis=1), 10, 7, 20000)
    assert grid[np.argmin(y)] >= 14000                           # a plain argmin of the whole run lands on the late dip: the guard is what prevents that
    ep_clean, _ = s.cell_epoch(LOG, np.stack([u_curve(rng, False, False) for _ in range(6)], axis=1), cols, A)
    assert abs(ep_clean - ep) <= 40                              # with and without the late trouble the pick is the same


def test_no_overfitting_found_is_reported():
    mat = np.stack([0.3 + np.exp(-LOG / 5000)] * 3, axis=1)     # still falling at the end of the window
    ep, found = s.cell_epoch(LOG, mat, np.arange(3), A)
    assert not found and ep >= 2900                              # the end of the search window, flagged


def test_snapping_goes_to_the_nearest_checkpoint():
    grid = list(range(50, 1001, 50))
    assert s.snap(239, grid) == 250 and s.snap(224, grid) == 200 and s.snap(225, grid) == 200 and s.snap(1e5, grid) == 1000


def scores(rng, means, sd, n_it=6):
    idx = pd.MultiIndex.from_product([["E1", "E2"], ["L1", "L2", "L3", "L4"], range(1, n_it + 1)], names=["embedding", "label", "iteration"])
    return pd.DataFrame({l: means.get(l, 0.8) + rng.normal(0, sd, len(idx)) for l in s.LOSSES}, index=idx)


def test_loss_scopes():
    rng = np.random.default_rng(1)
    sc = scores(rng, {"weighted_bce": 0.82}, 0.01)               # weighted_bce is the better loss everywhere ...
    sc.loc[("E1", "L1"), "focal"] = 0.95 + rng.normal(0, 0.005, 6)  # ... except one cell where focal is far better
    g = s.choose(sc, "global", 2.0)
    assert set(g.values()) == {"weighted_bce"}
    h = s.choose(sc, "hierarchical", 2.0)
    assert h[("E1", "L1")] == "focal" and sum(v != "weighted_bce" for v in h.values()) == 1   # the global loss everywhere, the one real exception only
    c = s.choose(sc, "cell", 2.0)
    assert c[("E1", "L1")] == "focal" and set(c[k] for k in c if k != ("E1", "L1")) <= set(s.LOSSES)
    by_label = s.choose(sc, "label", 2.0)
    assert by_label[("E1", "L2")] == by_label[("E2", "L2")]     # one loss per label, whatever the embedding
    tie = pd.DataFrame({l: 0.5 for l in s.LOSSES}, index=sc.index)
    assert set(s.choose(tie, "cell", 2.0).values()) == {"bce"}   # an exact tie goes to the first loss


def test_rounding_up_to_the_next_hundred():
    grid = list(range(50, 1001, 50))
    up = lambda e: s.round_epoch(e, "up", 100, grid)
    assert up(153) == 200 and up(142.5) == 200 and up(130) == 200 and up(180) == 200       # everything in (100, 200] goes to 200
    assert up(200) == 200 and up(100) == 100 and up(200.01) == 300 and up(199.9999999) == 200 and up(201) == 300  # an exact multiple stays; just above it goes up (the step function)
    assert up(1e5) == 1000 and up(10) == 100                                                # beyond the last checkpoint: the last one; below the first multiple: 100
    assert s.round_epoch(153, "nearest", 100, grid) == 150 and s.round_epoch(153, "up", 0, grid) == 150   # 'nearest' and unit 0 keep the old behaviour


def test_epoch_per_embedding_and_losses_read_at_that_epoch():
    grid = list(range(50, 1001, 50))
    pc = pd.DataFrame({"embedding": ["A"] * 4 + ["B"] * 4, "loss_type": s.LOSSES * 2, "epoch": [100, 120, 140, 160, 300, 320, 340, 360]})
    raw, rounded = s.epochs_by_embedding(pc, "mean", grid, "embedding")
    assert raw == {"A": 130.0, "B": 330.0} and rounded == {"A": 200.0, "B": 400.0}       # each embedding: the average of ITS OWN four losses, rounded UP to the next 100
    _, snapped = s.epochs_by_embedding(pc, "mean", grid, "embedding", rounding="nearest")
    assert snapped == {"A": 150.0, "B": 350.0}                                           # the earlier rule: the nearest checkpoint
    _, glob = s.epochs_by_embedding(pc, "mean", grid, "global")
    assert glob == {"A": 300.0, "B": 300.0}                                              # global scope: the average of all eight (230 -> up), one epoch for both
    _, hand = s.epochs_by_embedding(pc, "mean", grid, "embedding", manual=230)
    assert hand == {"A": 250.0, "B": 250.0}                                              # a manual epoch overrides the scope and is NOT rounded up (nearest checkpoint)
    # the MCC of the losses is read at each embedding's own epoch: A's best loss at 150 is bce, B's best loss at 350 is focal (and the opposite at the other epoch)
    rows = []
    for emb, best_at in (("A", {150: "bce", 350: "focal"}), ("B", {150: "bce", 350: "focal"})):
        for it in (1, 2, 3):
            for ep in (150, 350):
                for loss in ("bce", "focal"):
                    rows.append({"embedding": emb, "loss_type": loss, "iteration": it, "epoch": ep, "label": "L", "mcc": 0.9 if loss == best_at[ep] else 0.7})
    ck = pd.DataFrame(rows)
    sc = s.scores_per_embedding(ck, snapped, 0)
    m = s.choose(sc, "cell", 2.0)
    assert m[("A", "L")] == "bce" and m[("B", "L")] == "focal"                          # A read at 150, B read at 350
    sc_all150 = s.scores_per_embedding(ck, {"A": 150.0, "B": 150.0}, 0)
    assert s.choose(sc_all150, "cell", 2.0)[("B", "L")] == "bce"                         # at one common epoch B would have been judged at the wrong place


def test_selecting_among_equal_losses_is_optimistic():
    rng = np.random.default_rng(2)
    sc = scores(rng, {}, 0.05)                                   # four identical losses: any difference is noise
    held, fixed, avg = s.leave_one_split_out(sc, "cell", 2.0)
    insample = s.score_map(sc, s.choose(sc, "cell", 2.0))
    assert insample - held > 0.01                                # picking the best of four noisy cells flatters the in-sample score
    assert abs(held - avg) < 0.03                                # and on held-out splits the pick is no better than average


if __name__ == "__main__":
    for t in (test_guard_ignores_late_dips_and_spikes, test_no_overfitting_found_is_reported, test_snapping_goes_to_the_nearest_checkpoint, test_rounding_up_to_the_next_hundred,
              test_loss_scopes, test_epoch_per_embedding_and_losses_read_at_that_epoch, test_selecting_among_equal_losses_is_optimistic):
        t()
        print("ok", t.__name__)
