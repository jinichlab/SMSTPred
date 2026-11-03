"""sam_vs_up and its report: the Nadeau-Bengio test inflates the variance as it should, the vector loader keeps only what the H5 holds, the arms' metrics use the fixed 0.5 threshold;
classical_ml: the seven classifiers of Fig S4 fit and predict on a toy problem. Offline, seconds."""
import tempfile
from pathlib import Path

import h5py
import numpy as np

import bakeoff as b
import classical_ml as c
import sam_vs_up as s
import sam_vs_up_report as r


def test_nadeau_bengio_formula_and_inflation():
    rng = np.random.default_rng(0)
    d = rng.normal(0.02, 0.03, 16)
    t, p = r.nadeau_bengio(d, 0.3)
    n = len(d)
    assert np.isclose(t, d.mean() / np.sqrt((1 / n + 0.3) * d.var(ddof=1)))              # the formula
    from scipy import stats
    assert p > stats.ttest_rel(d, np.zeros(n)).pvalue                                    # always more conservative than the naive paired t-test
    assert r.nadeau_bengio(d, 0.0)[1] == stats.ttest_1samp(d, 0).pvalue or abs(r.nadeau_bengio(d, 0.0)[1] - stats.ttest_1samp(d, 0).pvalue) < 1e-12   # rho = 0 gives the naive test
    assert r.nadeau_bengio(np.zeros(16), 0.3)[1] == 1.0                                  # no difference at all


def test_load_vectors_keeps_only_what_is_there():
    with tempfile.TemporaryDirectory() as tmp:
        f = Path(tmp) / "v.h5"
        with h5py.File(f, "w") as h:
            for e in ("A", "B"):
                h[e] = np.arange(4, dtype=np.float16)
        v = s.load_vectors(f, {"A", "B", "C"})
        assert set(v) == {"A", "B"} and v["A"].dtype == np.float32
        ok, X = s.stack(v, ["B", "C", "A"])
        assert ok == ["B", "A"] and X.shape == (2, 4)


def test_arm_metrics_use_the_fixed_threshold():
    rng = np.random.default_rng(1)
    y = (rng.random((300, 4)) < 0.3).astype(np.float32)
    p = np.clip(y * 0.4 + rng.normal(0.3, 0.1, y.shape), 0, 1).astype(np.float32)
    overall, per_label = {}, {}
    s.metrics_rows("X", y, p, 1, 1000, 100, 300, overall, per_label)
    assert np.isclose(overall["X"]["mcc"], b.all_metrics(y, p, 0.5)["mcc"]) and overall["X"]["n_train"] == 1000 and set(per_label["X"]) == set(b.LABELS)


def test_classical_classifiers_fit_and_predict():
    rng = np.random.default_rng(2)
    X = rng.normal(size=(240, 16)).astype(np.float32)
    y = ((X[:, :4] + rng.normal(0, 0.5, (240, 4))) > 0).astype(np.float32)
    tmpl = c.templates(2)
    assert list(tmpl) == c.CLASSIFIERS and len(tmpl) == 7
    for name in ("Logistic Regression", "SVM_Linear", "Random Forest"):
        p = c.fit_predict(name, tmpl[name], X[:180], y[:180], X[180:], 2)
        assert p.shape == (60, 4) and np.isfinite(p).all() and 0 <= p.min() and p.max() <= 1
        assert np.mean((p > 0.5) == (y[180:] > 0.5)) > 0.7, name                          # it learned the signal
    y1 = y.copy(); y1[:, 3] = 0
    p = c.fit_predict("SVM_RBF", tmpl["SVM_RBF"], X[:180], y1[:180], X[180:], 2)
    assert (p[:, 3] == 0).all()                                                           # a label without positives in train: probability 0, no crash


if __name__ == "__main__":
    for f in (test_nadeau_bengio_formula_and_inflation, test_load_vectors_keeps_only_what_is_there, test_arm_metrics_use_the_fixed_threshold, test_classical_classifiers_fit_and_predict):
        f()
        print("ok", f.__name__)
