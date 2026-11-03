#!/usr/bin/env python3
"""
Choose the training epoch of every (embedding, label) from the validation loss alone: plain BCE, no MCC, no loss selection.

  python select_epoch_per_label.py [--dir train_val_loss_output] [--set sam_full] [--loss bce] [--splits 6] [--out epoch_selection_per_label_output]

The rule (the overfitting guard of select_epoch_and_loss.py, applied to one label's curve):
  For each (embedding, label) the BCE validation loss of that label (jointly trained 4-label MLP, as recorded by train_val_loss_curves.py), averaged over the
  cluster-disjoint splits and smoothed, is followed from the first epoch on until OVERFITTING SETS IN: the first moment it has stayed more than --tol (3 %) above its
  running minimum for --patience (100) epochs. The epoch of the minimum before that moment is the epoch of that (embedding, label). Later epochs are never looked at
  (late spikes and dips cannot be picked); --max-search (3,000) is a second hard limit. Epochs are not moved to a checkpoint: nothing here needs the MCC.

In the bake-off the embedding's MLP is still ONE network with a shared trunk and four outputs, trained once with BCE: label j's output is read from the network as it was
at the epoch chosen for label j (per-task early stopping; the network is trained to the largest of the four epochs and snapshotted on the way), and label j's
threshold is tuned on that snapshot's training predictions. That is exactly the state in which the sweep measured the label's validation loss.

Output (<out>): selection_per_label_<set>.json ({"epoch_by_embedding_label": {embedding: {label: epoch}}, "loss_map": {embedding: {label: "bce"}}, intervals, ...}),
epoch_per_label_<set>.csv, report, and epoch_selection_per_label_<set>.png (train/validation loss of every embedding x label with the chosen epoch).
The bake-off (mmseqs_clustering/sam_classifier_comparison.py, LOSS_MAP_JSON=selection_per_label_<set>.json) reads it directly.
"""
import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

import select_epoch_and_loss as pick   # smooth_curve and onset_epoch: the same guard as the loss-averaged picker

HERE = Path(__file__).resolve().parent
LABELS = ["Metabolite", "Protein", "RNA", "DNA"]


def label_cells(curves, loss):
    """{(embedding, label): (epochs, matrix [epochs x splits] of the validation loss)} of one loss."""
    cells = {}
    for (emb, lab), g in curves[curves.loss_type == loss].groupby(["embedding", "label"]):
        piv = g.pivot_table(index="epoch", columns="iteration", values="val_loss")
        cells[(emb, lab)] = (piv.index.values.astype(float), piv.values)
    return cells


def cell_epoch(epochs, mat, cols, a):
    """(epoch before overfitting sets in, whether an onset was found) of the split-averaged validation loss of one cell, using the columns `cols` (the splits)."""
    grid, y = pick.smooth_curve(epochs, mat[:, cols].mean(axis=1), a.step, a.smooth, a.max_search)
    return pick.onset_epoch(grid, y, a.tol, max(1, int(round(a.patience / a.step))))


def select(cells, a, n_boot, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for (emb, lab), (epochs, mat) in cells.items():
        n_it = mat.shape[1]
        ep, found = cell_epoch(epochs, mat, np.arange(n_it), a)
        boot = np.array([cell_epoch(epochs, mat, rng.integers(0, n_it, n_it), a)[0] for _ in range(n_boot)]) if n_boot else np.array([ep])   # --bootstrap 0: no interval (it collapses to the epoch)
        rows.append({"embedding": emb, "label": lab, "epoch": ep, "onset_found": found, "ci95_lo": float(np.percentile(boot, 2.5)), "ci95_hi": float(np.percentile(boot, 97.5)),
                     "n_splits": n_it})
    return pd.DataFrame(rows)


def plot(cells, cells_train, table, loss, out, a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    embs = list(dict.fromkeys(table.embedding))
    labs = [l for l in LABELS if l in set(table.label)]
    fig, ax = plt.subplots(len(embs), len(labs), figsize=(4.2 * len(labs), 2.7 * len(embs)), sharex=True, squeeze=False)
    for i, e in enumerate(embs):
        for j, l in enumerate(labs):
            epochs, va = cells[(e, l)]
            tr = cells_train[(e, l)][1]
            m = epochs <= a.max_search * 1.5
            for k in range(va.shape[1]):
                ax[i, j].plot(epochs[m], va[m, k], color="#d62728", alpha=0.18, lw=0.7)
            ax[i, j].plot(epochs[m], va[m].mean(axis=1), color="#d62728", lw=1.5, label="validation (mean of splits; thin = each split)")
            ax[i, j].plot(epochs[m], tr[m].mean(axis=1), color="#1f77b4", lw=1.2, label="train (mean of splits)")
            r = table[(table.embedding == e) & (table.label == l)].iloc[0]
            ax[i, j].axvspan(r.ci95_lo, r.ci95_hi, color="k", alpha=0.08)
            ax[i, j].axvline(r.epoch, color="k", ls="--", lw=1.1)
            ax[i, j].set_xscale("log"); ax[i, j].set_yscale("log")
            ax[i, j].set_title(f"{e} | {l}: epoch {r.epoch:.0f}", fontsize=9)
            ax[i, j].set_xlim(1, a.max_search)
            if j == 0:
                ax[i, j].set_ylabel(f"{loss} loss")
            if i == len(embs) - 1:
                ax[i, j].set_xlabel("epoch")
    ax[0, 0].legend(fontsize=6)
    fig.suptitle(f"{loss.upper()} train/validation loss per label and embedding, cluster-disjoint splits; dashed = chosen epoch (minimum before overfitting), band = 95 % bootstrap interval over splits", fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)


def run_loss(a, curves, loss, out):
    """The whole selection for one loss: report, json, csv, figure into `out`. -> (table, {embedding: {label: epoch}}, {embedding: {label: [lo, hi]}}, n_splits)."""
    cells = label_cells(curves[curves.loss_type == loss], loss)
    cells_train = {}
    for (emb, lab), g in curves[curves.loss_type == loss].groupby(["embedding", "label"]):
        piv = g.pivot_table(index="epoch", columns="iteration", values="train_loss")
        cells_train[(emb, lab)] = (piv.index.values.astype(float), piv.values)
    n_it = next(iter(cells.values()))[1].shape[1]
    table = select(cells, a, a.bootstrap)

    lines = []
    say = lambda s="": (print(s, flush=True), lines.append(s))
    say(f"{a.set}: loss {loss}, {len(table.embedding.unique())} embeddings x {len(table.label.unique())} labels, {n_it} splits" + (f" (the first {a.splits} of the sweep; the rest are left out of the selection)" if a.splits else ""))
    say(f"epoch = minimum of the split-averaged, smoothed validation loss of the label before overfitting sets in (> {a.tol:.0%} above the running minimum for {a.patience} epochs; search limited to {a.max_search} epochs)")
    wide = table.pivot(index="embedding", columns="label", values="epoch")[[l for l in LABELS if l in set(table.label)]]
    ci = table.assign(txt=lambda t: t.ci95_lo.round().astype(int).astype(str) + "-" + t.ci95_hi.round().astype(int).astype(str)).pivot(index="embedding", columns="label", values="txt")
    say("\nchosen epoch per embedding x label (95 % bootstrap interval over the splits):")
    for emb in wide.index:
        say(f"   {emb:16s} " + "  ".join(f"{l} {wide.loc[emb, l]:5.0f} [{ci.loc[emb, l]}]" for l in wide.columns))
    say(f"\nno overfitting onset found in {int((~table.onset_found).sum())} of {len(table)} cells (their epoch is the minimum of the search window)")
    say("mean over the four labels of each embedding (for comparison with the earlier per-embedding rule): " + ", ".join(f"{e} {wide.loc[e].mean():.0f}" for e in wide.index))
    say("by label (mean over embeddings): " + ", ".join(f"{l} {wide[l].mean():.0f}" for l in wide.columns))
    noisy = table[(table.ci95_hi - table.ci95_lo) > 1.0 * table.epoch]
    if len(noisy):
        say("NOTE: the interval is wider than the epoch itself in " + ", ".join(f"{r.embedding}/{r.label}" for r in noisy.itertuples()) + " (a rare label has few validation positives per split: its minimum is poorly determined)")

    out.mkdir(parents=True, exist_ok=True)
    by = {e: {l: float(wide.loc[e, l]) for l in wide.columns} for e in wide.index}
    iv = {e: {l: [float(table[(table.embedding == e) & (table.label == l)].ci95_lo.iloc[0]), float(table[(table.embedding == e) & (table.label == l)].ci95_hi.iloc[0])] for l in wide.columns} for e in wide.index}
    (out / f"selection_per_label_{a.set}.json").write_text(json.dumps({
        "set": a.set, "rule": "per-(embedding,label) epoch from the validation loss of one loss; no MCC", "loss": loss,
        "epoch_by_embedding_label": by, "epoch_interval_95_by_embedding_label": iv,
        "loss_map": {e: {l: loss for l in wide.columns} for e in wide.index},
        "selected_on_first_n_splits": a.splits or None, "n_splits_used": int(n_it), "overfitting_tolerance": a.tol, "patience": a.patience, "max_search": a.max_search,
    }, indent=1))
    table.to_csv(out / f"epoch_per_label_{a.set}.csv", index=False)
    (out / f"selection_report_per_label_{a.set}.txt").write_text("\n".join(lines) + "\n")
    if not a.no_plot:
        plot(cells, cells_train, table, loss, out / f"epoch_selection_per_label_{a.set}.png", a)
    say(f"wrote {out}/")
    return table, by, iv, int(n_it)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=str(HERE / "train_val_loss_output"))
    ap.add_argument("--set", default="sam_full")
    ap.add_argument("--loss", default="bce", help="the one loss used for every embedding and label; 'all' = each of the four losses in turn (<out>/<loss>/ and one combined json)")
    ap.add_argument("--out", default=str(HERE / "epoch_selection_per_label_output"))
    ap.add_argument("--splits", type=int, default=0, help="select on the first N splits only (0 = all); the bake-off keeps splits N+1.. as an unbiased check")
    ap.add_argument("--tol", type=float, default=0.03)
    ap.add_argument("--patience", type=int, default=100)
    ap.add_argument("--max-search", type=int, default=3000)
    ap.add_argument("--smooth", type=int, default=7)
    ap.add_argument("--step", type=int, default=10)
    ap.add_argument("--bootstrap", type=int, default=200)
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args()

    d = Path(a.dir)
    files = sorted(d.glob(f"{a.set}__*__curves.csv.gz"))
    if not files:
        sys.exit(f"no {a.set}__*__curves.csv.gz in {d}: run train_val_loss_curves.py first")
    cols = ["embedding", "loss_type", "iteration", "epoch", "label", "train_loss", "val_loss"]
    curves = pd.concat([pd.read_csv(f, usecols=cols) for f in files], ignore_index=True)
    if a.splits:
        curves = curves[curves.iteration <= a.splits]
    losses = pick.LOSSES if a.loss == "all" else [a.loss]
    missing = [l for l in losses if not (curves.loss_type == l).any()]
    if missing:
        sys.exit(f"no rows of loss {missing} in {d}")
    out = Path(a.out)
    results = {}
    for loss in losses:
        results[loss] = run_loss(a, curves, loss, out / loss if a.loss == "all" else out)
    if a.loss == "all":
        # one file for the bake-off (mmseqs_clustering/sam_classifier_comparison.py): every loss trained as its own network, each label read at ITS epoch of THAT loss
        (out / f"selection_per_label_all_losses_{a.set}.json").write_text(json.dumps({
            "set": a.set, "rule": "per-(loss,embedding,label) epoch from that loss's validation loss; the loss per label is chosen AFTER the bake-off, out of sample (assemble_best_loss.py)",
            "losses": losses, "epoch_by_loss_embedding_label": {l: results[l][1] for l in losses},
            "epoch_interval_95_by_loss_embedding_label": {l: results[l][2] for l in losses},
            "selected_on_first_n_splits": a.splits or None, "n_splits_used": results[losses[0]][3],
            "overfitting_tolerance": a.tol, "patience": a.patience, "max_search": a.max_search}, indent=1))
        print(f"\nwrote {out}/selection_per_label_all_losses_{a.set}.json", flush=True)


if __name__ == "__main__":
    main()
