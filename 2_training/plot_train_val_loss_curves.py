#!/usr/bin/env python3
"""
Train vs. validation LOSS vs epoch grids -- overfitting counterparts of the
MCC-vs-epoch Figures S2 (SAM full) and S5 (UP ProstT5); the grids are Figures S3 and S6.

Reads the per-task curves written by train_val_loss_curves.py and writes, into
train_val_loss_output/:
  full_train_val_loss_grid.png       SAM full      rows = embedding, cols = loss type
  balanced_train_val_loss_grid.png   SAM balanced  (only if balanced curves exist; not in the manuscript)
  UP_prostt5_train_val_loss_grid.png UP ProstT5    single row,       cols = loss type
  per_label/<set>_<loss>_per_label_grid.png   rows = embedding, cols = label (one
                                               file per loss; the label-resolved view)
  train_val_loss_summary.csv         epoch of minimum validation loss, loss there,
                                     loss at the end, loss at the chosen epoch
  train_val_loss_summary_by_iteration.csv   same minimum, per split (val sets range 36-5,771
                                     proteins, so the mean curve can hide big split-to-split spread)

Each panel: mean over the 6 cluster-disjoint iterations as the line, the band spanning
the min to max over the splits (a +/- SD band would run below zero on the log axis because
the validation sets differ so much in size), plus one faint line per split for validation;
train solid, validation dashed; log x and log y (training loss collapses toward 0 while
validation loss rises, so linear y hides one or the other). Dot = epoch of
minimum mean validation loss. Grey dotted vertical = the epoch used in the
benchmark (300 full set, 100 balanced set, 400 for the UP-trained ProstT5). Loss
panels are the optimized objective, i.e. summed over the four labels; each
loss type keeps its own y-scale (weighted variants are not comparable to bce).
Pure CSV -> PNG, no GPU.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatter
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
IN_DIR = HERE / "train_val_loss_output"     # overridable with --dir
LABELS = ["DNA", "RNA", "Protein", "Metabolite"]               # S2/S5 column order
LOSSES = ["bce", "weighted_bce", "sqrt_weighted_bce", "focal"]
EMBS = ["ProstT5", "SaProt", "ProtT5", "ESM-IF", "ProtT5+ESM-IF", "ProstT5+ESM-IF"]
C_TRAIN, C_VAL = "#4578aa", "#c14445"


class PlainLogFormatter(LogFormatter):
    """Log-axis tick labels as plain numbers (0.4, 1, 2) instead of 4e-01 / 4x10^-1.
    The parent still decides WHICH ticks get a label (minor ones only on narrow spans)."""
    def _num_to_string(self, x, vmin, vmax):
        return f"{x:g}"

# Epoch used in the manuscript's benchmark per set: every full-set MLP at 300 epochs (Fig 2A / Table 1), every
# balanced-set MLP at 100 (supplementary, sam_full_vs_balanced_classifier_comparison.py presets), and the
# UP-trained ProstT5 at 400 (Fig 2B).
CHOSEN = {
    "full":     {e: 300 for e in ("ProtT5", "ProstT5", "ESM-IF", "ProtT5+ESM-IF", "ProstT5+ESM-IF", "SaProt")},
    "balanced": {e: 100 for e in ("ProtT5", "ProstT5", "ESM-IF", "ProtT5+ESM-IF", "ProstT5+ESM-IF", "SaProt")},
    "UP":       {"ProstT5": 400},
}


def load(prefix):
    files = sorted(IN_DIR.glob(f"{prefix}__*__curves.csv.gz"))
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True) if files else None


def per_iter(df, by_label):
    """Per-iteration curves. by_label=False sums the four labels first (= the objective)."""
    keys = ["embedding", "loss_type", "iteration", "epoch"] + (["label"] if by_label else [])
    return df.groupby(keys)[["train_loss", "val_loss"]].sum().reset_index()


def mean_sd(d, by_label):
    keys = ["embedding", "loss_type", "epoch"] + (["label"] if by_label else [])
    g = d.groupby(keys)[["train_loss", "val_loss"]].agg(["mean", "std"])
    g.columns = [f"{a}_{b}" for a, b in g.columns]
    return g.reset_index()


def draw(ax, g, raw, chosen):
    x = g["epoch"].values
    for _, ri in raw.groupby("iteration"):          # one faint line per split: val sets differ a lot in size
        ax.plot(ri["epoch"].values, ri["val_loss"].values, color=C_VAL, lw=0.5, alpha=0.25, zorder=2)
    for col, color, ls, lab in (("train_loss", C_TRAIN, "-", "train"),
                                ("val_loss", C_VAL, "--", "validation")):
        m = g[f"{col}_mean"].values
        rng = raw.groupby("epoch")[col].agg(["min", "max"]).reindex(x)   # spread over the 6 splits (log-axis safe)
        ax.plot(x, m, color=color, ls=ls, lw=1.6, label=lab, zorder=3)
        ax.fill_between(x, rng["min"].values, rng["max"].values, color=color, alpha=0.12, lw=0, zorder=1)
    vm = g["val_loss_mean"].values
    pk = int(vm.argmin())
    ax.plot(x[pk], vm[pk], "o", color=C_VAL, ms=5, zorder=4)
    if chosen is not None:
        ax.axvline(chosen, color="0.45", ls=":", lw=1.2, zorder=2)
    ax.set_xscale("log"); ax.set_yscale("log")
    for which in ("major", "minor"):                # plain numbers (0.4, 1, 2) instead of 4x10^-1
        getattr(ax.yaxis, f"set_{which}_formatter")(PlainLogFormatter(labelOnlyBase=False, minor_thresholds=(2, 0.5)))
    ax.grid(alpha=0.3, which="both")


def grid(df, rows, cols, col_field, by_label, outpath, set_key, sharey):
    raw_all = per_iter(df, bool(by_label) or col_field == "label")
    g_all = mean_sd(raw_all, bool(by_label) or col_field == "label")
    rows = [r for r in rows if r in set(g_all["embedding"])]
    nr, nc = len(rows), len(cols)
    fig, axes = plt.subplots(nr, nc, figsize=(4.4 * nc, 3.0 * nr + 1.0), squeeze=False,
                             sharex=True, sharey=sharey)
    for r, emb in enumerate(rows):
        for c, colv in enumerate(cols):
            ax = axes[r][c]
            m = lambda t: (t.embedding == emb) & (t[col_field] == colv)
            sub, raw = g_all[m(g_all)], raw_all[m(raw_all)]
            if by_label:
                sub, raw = sub[sub.loss_type == by_label], raw[raw.loss_type == by_label]
            if sub.empty:
                ax.set_visible(False); continue
            draw(ax, sub.sort_values("epoch"), raw.sort_values("epoch"), CHOSEN[set_key].get(emb))
            if r == 0:
                ax.set_title(colv, fontsize=12, fontweight="bold")
            if c == 0:
                ax.set_ylabel(f"{emb}\nloss", fontsize=9)
            if r == nr - 1:
                ax.set_xlabel("epoch (log)")
    n_iters = df["iteration"].nunique()
    axes[0][-1].legend(fontsize=7, loc="upper right", title=f"(n={n_iters})")
    fig.tight_layout()
    fig.savefig(outpath, dpi=140)
    plt.close(fig)
    print("wrote", outpath)


def summarize(df, set_key):
    raw = per_iter(df, by_label=False)
    g = mean_sd(raw, by_label=False)
    rows, it_rows = [], []
    for (emb, loss), s in g.groupby(["embedding", "loss_type"]):
        s = s.sort_values("epoch")
        pk = int(s["val_loss_mean"].values.argmin())
        ch = CHOSEN[set_key].get(emb)
        at = s[s.epoch == ch] if ch is not None else s.iloc[0:0]
        last = s.iloc[-1]
        rows.append({"dataset": set_key, "embedding": emb, "loss_type": loss,
                     "epoch_min_val": int(s["epoch"].values[pk]),
                     "val_min": s["val_loss_mean"].values[pk],
                     "chosen_epoch": ch,
                     "val_at_chosen": at["val_loss_mean"].iloc[0] if len(at) else np.nan,
                     "train_at_chosen": at["train_loss_mean"].iloc[0] if len(at) else np.nan,
                     "last_epoch": int(last["epoch"]),
                     "val_last": last["val_loss_mean"], "train_last": last["train_loss_mean"],
                     "val_last_over_min": last["val_loss_mean"] / s["val_loss_mean"].values[pk]})
    for (emb, loss, it), s in raw.groupby(["embedding", "loss_type", "iteration"]):
        s = s.sort_values("epoch")
        pk = int(s["val_loss"].values.argmin())
        it_rows.append({"dataset": set_key, "embedding": emb, "loss_type": loss, "iteration": it,
                        "epoch_min_val": int(s["epoch"].values[pk]), "val_min": s["val_loss"].values[pk],
                        "val_last": s["val_loss"].values[-1], "train_last": s["train_loss"].values[-1]})
    return rows, it_rows


def main():
    global IN_DIR
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=None, help="curves directory (default train_val_loss_output)")
    a = ap.parse_args()
    if a.dir:
        IN_DIR = Path(a.dir).resolve()
    (IN_DIR / "per_label").mkdir(parents=True, exist_ok=True)
    summary, by_iter = [], []
    for key, prefix, name in (("full", "sam_full", "full_train_val_loss_grid.png"),
                              ("balanced", "sam_balanced", "balanced_train_val_loss_grid.png"),
                              ("UP", "up", "UP_prostt5_train_val_loss_grid.png")):
        df = load(prefix)
        if df is None:
            print(f"[skip] no curves for {prefix}"); continue
        grid(df, EMBS, LOSSES, "loss_type", False, IN_DIR / name, key, sharey="col")
        for loss in LOSSES:
            if (df.loss_type == loss).any():
                grid(df, EMBS, LABELS, "label", loss, IN_DIR / "per_label" / f"{prefix}_{loss}_per_label_grid.png",
                     key, sharey="row")
        s_rows, i_rows = summarize(df, key)
        summary += s_rows; by_iter += i_rows
    pd.DataFrame(summary).to_csv(IN_DIR / "train_val_loss_summary.csv", index=False)
    pd.DataFrame(by_iter).to_csv(IN_DIR / "train_val_loss_summary_by_iteration.csv", index=False)
    print("wrote", IN_DIR / "train_val_loss_summary.csv", "and _by_iteration.csv")


if __name__ == "__main__":
    main()
