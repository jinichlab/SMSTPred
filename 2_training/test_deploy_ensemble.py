"""deploy_ensemble: whole clusters in exactly one fold, balanced by size and by label, deterministic; the applicability-domain cutoff, the calibration table and the ensemble mean / spread behave;
the saved members load into Stage 3's own model class and give the same output. Offline, seconds."""
import importlib.util
import tempfile
from pathlib import Path

import numpy as np
import torch

import bakeoff as b
import deploy_ensemble as d


def toy(n_clusters=300, seed=0):
    rng = np.random.default_rng(seed)
    sizes = np.minimum((rng.pareto(1.1, n_clusters) * 2 + 1).astype(int), 120)          # heavy-tailed cluster sizes
    cl = np.repeat(np.arange(n_clusters), sizes)
    y = np.zeros((len(cl), 4), np.float32)
    y[:, 0] = rng.random(len(cl)) < 0.2; y[:, 1] = rng.random(len(cl)) < 0.15; y[:, 2] = rng.random(len(cl)) < 0.6; y[:, 3] = rng.random(len(cl)) < 0.02
    y[:, 3][::40] = 1
    return cl, y


def test_pack_folds_whole_clusters_balanced_and_deterministic():
    cl, y = toy()
    fold, load, worst = d.pack_folds(cl, y, 6, 42)
    assert set(fold) == set(range(6))
    for c in np.unique(cl):
        assert len(set(fold[cl == c])) == 1                                              # a cluster is never split
    sizes = np.bincount(fold)
    assert sizes.max() / sizes.min() < 1.25, sizes                                       # balanced by protein count
    assert (np.array([y[fold == j].sum(0) for j in range(6)]) > 0).all()                 # every label occurs in every fold
    assert worst < 0.5
    assert np.array_equal(fold, d.pack_folds(cl, y, 6, 42)[0]) and not np.array_equal(fold, d.pack_folds(cl, y, 6, 43)[0])


def test_every_protein_is_validated_exactly_once():
    cl, y = toy(60)
    fold, _, _ = d.pack_folds(cl, y, 6, 1)
    seen = np.zeros(len(cl), int)
    for j in range(6):
        va, tr = fold == j, fold != j
        assert not (va & tr).any() and (va | tr).all()
        seen += va
    assert (seen == 1).all()


def test_audit_rotation_is_honest():
    cl, y = toy(60)
    fold, _, _ = d.pack_folds(cl, y, 6, 1)
    tested = np.zeros(len(cl), int)
    for j in range(6):
        te, va = fold == j, fold == (j + 1) % 6
        tr = ~(te | va)
        assert not (te & va).any() and not (te & tr).any() and not (va & tr).any()           # the test fold is neither trained on nor used for the epoch
        tested += te
    assert (tested == 1).all()                                                                # every protein is tested exactly once


def test_ad_cutoff_calibration_and_ensemble():
    X = np.vstack([np.eye(8, dtype=np.float32)] * 2 + [np.random.default_rng(0).normal(size=(4, 8)).astype(np.float32)])
    cut, nn = d.ad_cutoff_from(X, 1)
    assert np.allclose(nn[:16], 1.0) and cut <= 1.0 and nn.shape == (20,)                 # a duplicated vector has nearest neighbour similarity 1 (never itself)
    y = (np.random.default_rng(1).random((200, 4)) < 0.3).astype(np.float32)
    cal = d.calibration_table(y, np.clip(y * 0.6 + 0.2, 0, 1))
    assert set(cal.label) == set(b.LABELS) and cal.n.sum() == 800
    torch.manual_seed(0)
    models = [b.MultilabelNN(8, 4).to(b.DEVICE).eval() for _ in range(3)]
    mean, sd = d.ensemble_predict(models, X)
    assert mean.shape == sd.shape == (20, 4) and (sd > 0).all() and ((mean >= 0) & (mean <= 1)).all()


def test_saved_members_load_into_stage3_model_class():
    spec = importlib.util.spec_from_file_location("predict_stage3", Path(__file__).resolve().parents[1] / "3_predictions" / "predict.py")
    stage3 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(stage3)
    m = b.MultilabelNN(16, 4)
    with tempfile.TemporaryDirectory() as tmp:
        torch.save({k: v.cpu() for k, v in m.state_dict().items()}, Path(tmp) / "m.pt")
        s3 = stage3.MultilabelNN(16, 4, [512, 256])
        s3.load_state_dict(torch.load(Path(tmp) / "m.pt", map_location="cpu"))
    x = torch.randn(5, 16)
    assert torch.allclose(m.eval()(x), s3.eval()(x))                                      # same layers, same output: a drop-in for Stage 3


if __name__ == "__main__":
    for f in (test_pack_folds_whole_clusters_balanced_and_deterministic, test_every_protein_is_validated_exactly_once, test_audit_rotation_is_honest, test_ad_cutoff_calibration_and_ensemble, test_saved_members_load_into_stage3_model_class):
        f()
        print("ok", f.__name__)
