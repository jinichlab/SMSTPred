"""bakeoff.py: early stopping keeps the best epoch's weights while training goes on, the primary threshold is a fixed 0.5 (HMM: validation-tuned), mcc_valtuned is recorded next to it,
the data loader reproduces the benchmark set. Offline except the last test (it needs a release folder: $RELEASE_DIR)."""
import os

import numpy as np

import bakeoff as b


def test_early_stopper():
    s = b.EarlyStopper(patience=3)
    kept = [s.update(i + 1, v) for i, v in enumerate([5, 4, 3, 3.5, 3.4, 3.6, 2.0, 1.0])]      # best 3 at epoch 3; fires at epoch 6; the later 2.0 / 1.0 can no longer be kept
    assert kept == [True, True, True, False, False, False, False, False] and s.best_epoch == 3 and s.stop_epoch == 6
    s2 = b.EarlyStopper(patience=3)
    for i, v in enumerate([3, 2, 1, 0.5]):
        s2.update(i + 1, v)
    assert s2.best_epoch == 4 and s2.stop_epoch is None


def synthetic(n=300, d=24, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, d)).astype(np.float32)
    y = ((X @ rng.normal(size=(d, 4)) + rng.normal(0, 2.5, (n, 4))) > np.array([0, 1.0, -0.5, 2.0])).astype(np.float32)   # noisy, so the validation loss turns up
    return X, y


def test_kept_weights_are_the_best_epoch_and_all_epochs_are_trained():
    X, y = synthetic()
    Xtr, ytr, Xva, yva = X[:200], y[:200], X[200:], y[200:]
    model, hist, stop = b.train_early_stop(Xtr, ytr, Xva, yva, epochs=300, patience=15, lr=1e-3, seed=3)
    assert hist["train"].shape == (300, 4) and np.isfinite(hist["val"]).all()                  # ALL epochs ran, the curves are complete
    assert stop.stop_epoch is not None and stop.stop_epoch < 300 and abs(hist["val"][stop.best_epoch - 1].sum() - stop.best) < 1e-5
    ref, _, _ = b.train_early_stop(Xtr, ytr, Xva, yva, epochs=stop.best_epoch, patience=10 ** 6, lr=1e-3, seed=3)          # the same seed trained to exactly the best epoch
    assert np.allclose(b.predict(model, Xva), b.predict(ref, Xva), atol=1e-6)


def test_primary_threshold_is_half_and_valtuned_is_recorded():
    rng = np.random.default_rng(1)
    y_va = (rng.random((400, 4)) < 0.3).astype(np.float32)
    s_va = np.clip(y_va * 0.6 + rng.normal(0.2, 0.1, y_va.shape), 0, 1).astype(np.float32)          # on validation positives ~0.8, negatives ~0.2
    y_te = (rng.random((400, 4)) < 0.3).astype(np.float32)
    s_te = np.clip(y_te * 0.5 + rng.normal(0.2, 0.1, y_te.shape), 0, 1).astype(np.float32)
    overall, per_label, scores = {}, {}, {}
    b.score_method("M", s_va, s_te, y_va, y_te, overall, per_label, scores)
    assert np.isclose(overall["M"]["mcc"], b.all_metrics(y_te, s_te, 0.5)["mcc"])                    # primary = fixed 0.5
    thr = b.tune_thresholds(y_va, s_va)
    assert np.allclose(scores["M__thr_valtuned"], thr) and np.isclose(overall["M"]["mcc_valtuned"], b.all_metrics(y_te, s_te, thr)["mcc"])
    assert not np.allclose(thr, 0.5) and all("mcc_valtuned" in per_label["M"][l] for l in b.LABELS)


def test_hmm_uses_the_validation_threshold_as_its_primary():
    rng = np.random.default_rng(2)
    y_va = (rng.random((400, 4)) < 0.3).astype(np.float32)
    s_va = np.where(y_va == 1, rng.normal(60, 10, y_va.shape), rng.normal(20, 10, y_va.shape))        # bit scores: positives ~60, negatives ~20
    y_te = (rng.random((400, 4)) < 0.3).astype(np.float32)
    s_te = np.where(y_te == 1, rng.normal(60, 10, y_te.shape), rng.normal(20, 10, y_te.shape))
    overall, per_label, scores = {}, {}, {}
    b.score_method("H", s_va, s_te, y_va, y_te, overall, per_label, scores, hmm=True)
    assert overall["H"]["mcc"] > 0.7 and np.isclose(overall["H"]["mcc"], overall["H"]["mcc_valtuned"])      # a 0.5 on bit scores would be meaningless; the tuned threshold gives a sensible MCC
    assert (scores["H__thr_valtuned"] > 25).all() and (scores["H__thr_valtuned"] < 60).all()               # between the two clouds


def test_loader_reproduces_the_benchmark_set():
    if not os.environ.get("RELEASE_DIR"):
        print("   (skipped: set RELEASE_DIR to a release folder)")
        return
    df, mats = b.load_data(b.release_dir())
    assert len(df) == len(set(df.Entry)) and df.Sequence.notna().all() and all(m.shape[0] == len(df) for m in mats.values())
    assert set(b.LABELS) <= set(df.columns) and mats["prott5"].shape[1] == 1024 and mats["prostt5"].shape[1] == 1024 and mats["saprot"].shape[1] == 1280 and mats["esmif"].shape[1] == 512
    print(f"   loader: {len(df)} benchmark proteins")


if __name__ == "__main__":
    for f in (test_early_stopper, test_kept_weights_are_the_best_epoch_and_all_epochs_are_trained, test_primary_threshold_is_half_and_valtuned_is_recorded,
              test_hmm_uses_the_validation_threshold_as_its_primary, test_loader_reproduces_the_benchmark_set):
        f()
        print("ok", f.__name__)
