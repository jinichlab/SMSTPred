#!/usr/bin/env python3
"""
Choose the training epoch and the loss of the MLP classifiers automatically, from the sweep that train_val_loss_curves.py records.

  python select_epoch_and_loss.py [--dir train_val_loss_output] [--set sam_full] [--out epoch_loss_selection_output]

The rule
  1. EPOCH. For every (embedding, loss) the validation loss (mean over the four labels, then over the splits) is followed from the first epoch on, smoothed, until
     OVERFITTING SETS IN: the first moment it has been more than --tol (3 %) above its running minimum for --patience (100) epochs in a row. The epoch of the minimum before that
     moment is that cell's epoch. Later epochs are never looked at, so the late phase of training, where the loss jumps up and down and can dip below the true minimum by accident,
     cannot be picked (the same as early stopping with patience); --max-search (3,000) is a second, hard limit. The epoch of an EMBEDDING is the AVERAGE of the epochs of its own
     four losses (--aggregate mean, or median), ROUNDED UP to the next multiple of --round-to (100): 153 -> 200, 200 -> 200 (--round nearest: the nearest MCC checkpoint instead;
     train_val_loss_curves.py scores the MCC every 50 epochs). Up, because the epoch where the validation loss bottoms out is where the model starts to become over-confident,
     while the MCC keeps improving a little longer (on the refreshed data ProstT5 gains +0.018 MCC on the selection splits 1-6, better on 6 of 6, and +0.021 on the unseen splits
     7-16 from 150 to 200 epochs although its loss minima average 142). The rule is a step function (142 -> 200, 201 -> 300). --epoch-scope global instead averages over every
     cell and uses one epoch for all embeddings (the earlier behaviour); embeddings differ in how long they need, which is why the default is per embedding.
  2. LOSS. At each embedding's epoch the loss with the highest mean MCC over the splits is chosen for each (embedding, label) cell (--loss-scope cell). Ties go to the first of
     bce, weighted_bce, sqrt_weighted_bce, focal. All four losses of an embedding are compared at the same epoch (equal training budget).

What it also tells you, because the sweep data say a per-cell loss choice is easily within noise (see 2_training/README.md): the bootstrap interval of E* over the splits; for each cell
the winner's margin against the paired standard error over the splits and whether it holds one checkpoint earlier or later; and a leave-one-split-out check of the whole rule (choose on
the other splits, score on the held-out one), which is the honest score: the in-sample score of a selected loss is optimistic. --loss-scope hierarchical (one global loss, a cell only
changes it when its gain beats --min-gain-se paired standard errors), label (one loss per label) and global are the safer variants; --epoch fixes the epoch by hand.

Input  (<dir>): <set>__<embedding>[__<loss>]__curves.csv.gz and ...__ckpt_mcc.csv.gz of train_val_loss_curves.py (columns as written there)
Output (<out>): selection_<set>.json, loss_map_<set>.json ({embedding: {label: loss}}), selection_cells_<set>.csv, epoch_per_cell_<set>.csv, selection_report_<set>.txt, epoch_selection_<set>.png
Pure numpy / pandas (matplotlib for the figure); a minute on the full sweep.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
LOSSES = ["bce", "weighted_bce", "sqrt_weighted_bce", "focal"]
LINES = []


def say(msg=""):
    print(msg, flush=True)
    LINES.append(msg)


# ---------------------------------------------------------------- 1. the epoch
def build_cells(curves):
    """{(embedding, loss): (epochs, matrix [epochs x splits])} of the validation loss averaged over the labels."""
    cells = {}
    for (emb, loss), g in curves.groupby(["embedding", "loss_type"]):
        piv = g.groupby(["epoch", "iteration"]).val_loss.mean().unstack("iteration")
        cells[(emb, loss)] = (piv.index.values.astype(float), piv.values)
    return cells


def smooth_curve(epochs, vals, step, width, max_epoch):
    """The curve on a regular grid (the logging is every epoch to 1,000 and every 10 after), smoothed by a centred moving average."""
    grid = np.arange(step, max_epoch + 1, step, dtype=float)
    y = np.interp(grid, epochs, vals)
    return grid, pd.Series(y).rolling(width, center=True, min_periods=1).mean().values


def onset_epoch(grid, y, tol, patience_pts):
    """Epoch of the minimum BEFORE overfitting sets in. Overfitting has set in at the first point from which the curve stays more than tol above its running minimum for
    patience_pts grid points in a row. -> (epoch, whether such an onset was found)."""
    running_min = np.minimum.accumulate(y)
    bad = (y > running_min * (1 + tol)).astype(int)
    if len(bad) >= patience_pts:
        run = np.convolve(bad, np.ones(patience_pts, dtype=int), mode="valid")
        hits = np.flatnonzero(run == patience_pts)
        if hits.size:
            return float(grid[int(np.argmin(y[:hits[0] + 1]))]), True
    return float(grid[int(np.argmin(y))]), False


def cell_epoch(epochs, mat, cols, a, max_epoch=None):
    vals = mat[:, cols].mean(axis=1)
    grid, y = smooth_curve(epochs, vals, a.step, a.smooth, a.max_search if max_epoch is None else max_epoch)
    return onset_epoch(grid, y, a.tol, max(1, int(round(a.patience / a.step))))


def snap(e, grid):
    grid = np.asarray(sorted(grid), dtype=float)
    return float(grid[int(np.argmin(np.abs(grid - e)))])


def round_epoch(e, how, unit, grid):
    """An averaged epoch -> the MCC checkpoint that is used. how 'up' (the default rule): the next multiple of `unit` (100), so 153 -> 200 and 200 -> 200 (an exact multiple stays),
    then onto the checkpoint grid (a multiple of 50 is one; beyond the last checkpoint it is the last). how 'nearest' (or unit 0): the nearest checkpoint, as before.
    Why up: the epoch where the validation loss bottoms out is where the model starts to become over-confident, while the MCC (threshold tuned on the training predictions) keeps
    improving a little longer, so the loss minimum is an early estimate of the best epoch. It is a step function: 142 -> 200, but 201 -> 300."""
    if how == "up" and unit:
        e = float(np.ceil(e / unit - 1e-9) * unit)
    return snap(e, grid)


def aggregate(values, how):
    return float(np.mean(values) if how == "mean" else np.median(values))


# ---------------------------------------------------------------- 2. the loss
def epochs_by_embedding(per_cell, how, grid, scope="embedding", manual=None, rounding="up", unit=100):
    """The epoch of every embedding -> ({embedding: aggregate of its cell epochs}, {embedding: the MCC checkpoint it is rounded to, see round_epoch}).
    scope 'embedding': each embedding's own, the average over its (four) losses of the epoch where that loss's validation loss bottoms out;
    scope 'global': one epoch for all, the average over every cell; manual: the same epoch for all (an epoch given by hand is moved to the nearest checkpoint, never rounded up)."""
    embs = list(dict.fromkeys(per_cell.embedding))
    if manual is not None:
        raw = {e: float(manual) for e in embs}
        return raw, {e: snap(v, grid) for e, v in raw.items()}
    elif scope == "embedding":
        raw = {e: aggregate(per_cell[per_cell.embedding == e].epoch, how) for e in embs}
    else:
        raw = {e: aggregate(per_cell.epoch, how) for e in embs}
    return raw, {e: round_epoch(v, rounding, unit, grid) for e, v in raw.items()}


def scores_per_embedding(ck, e_by, window):
    """The MCC scores of scores_at, each embedding read at ITS epoch; the result has the same layout, so choose() and the checks work unchanged."""
    return pd.concat([scores_at(ck[ck.embedding == e], ep, window) for e, ep in e_by.items()])


def scores_at(ck, epoch, window):
    """MCC per (embedding, label, split) and loss at the checkpoint (or the mean of the checkpoints within +-window epochs)."""
    sub = ck[(ck.epoch >= epoch - window) & (ck.epoch <= epoch + window)]
    sc = sub.pivot_table(index=["embedding", "label", "iteration"], columns="loss_type", values="mcc", aggfunc="mean")
    return sc[[l for l in LOSSES if l in sc.columns]]


def choose(sc, scope, k):
    """{(embedding, label): loss} from the scores sc (index embedding, label, split; columns losses)."""
    cell = sc.groupby(level=[0, 1]).mean()
    if scope == "cell":
        return cell.idxmax(axis=1).to_dict()  # idxmax keeps the first of equal maxima: ties go to the earlier loss of LOSSES
    best_global = cell.mean().idxmax()
    if scope == "global":
        return {c: best_global for c in cell.index}
    if scope == "label":
        by_label = cell.groupby(level=1).mean().idxmax(axis=1)
        return {c: by_label[c[1]] for c in cell.index}
    out = {}  # hierarchical: the global loss unless the cell's own best is better by more than k paired standard errors
    for c in cell.index:
        best = cell.loc[c].idxmax()
        out[c] = best_global
        if best != best_global:
            diff = (sc.loc[c][best] - sc.loc[c][best_global]).values
            se = diff.std(ddof=1) / np.sqrt(len(diff)) if len(diff) > 1 else np.inf
            if diff.mean() > k * se:
                out[c] = best
    return out


def score_map(sc, mapping):
    return float(np.mean([sc.loc[c][l].mean() for c, l in mapping.items()]))


def leave_one_split_out(sc, scope, k):
    """Choose with the other splits, score on the held-out one: the score a loss map made this way would have on new splits."""
    its = sorted(sc.index.get_level_values(2).unique())
    sel, fixed, avg = [], [], []
    for it in its:
        train = sc[sc.index.get_level_values(2) != it]
        m, g = choose(train, scope, k), choose(train, "global", k)
        held = sc[sc.index.get_level_values(2) == it].droplevel(2)
        sel.append(np.mean([held.loc[c][l] for c, l in m.items()]))
        fixed.append(np.mean([held.loc[c][g[c]] for c in m]))
        avg.append(held.mean(axis=1).mean())
    return float(np.mean(sel)), float(np.mean(fixed)), float(np.mean(avg))


def cell_table(sc, mapping, ck, e_by, window, grid, a):
    rows = []
    cell = sc.groupby(level=[0, 1]).mean()
    for (emb, lab), row in cell.iterrows():
        best = row.idxmax()
        rest = row.drop(best)
        run = rest.idxmax()
        diff = (sc.loc[(emb, lab)][best] - sc.loc[(emb, lab)][run]).values
        se = diff.std(ddof=1) / np.sqrt(len(diff)) if len(diff) > 1 else np.nan
        near = {}
        for step in (-a.checkpoint_step, a.checkpoint_step):
            e2 = e_by[emb] + step
            if e2 in grid:
                near[step] = scores_at(ck[ck.embedding == emb], e2, window).groupby(level=[0, 1]).mean().loc[(emb, lab)].idxmax()
        rows.append({"embedding": emb, "label": lab, "chosen_loss": mapping[(emb, lab)], "best_loss": best, "best_mcc": row[best], "runner_up": run,
                     "margin": row[best] - row[run], "paired_se": se, "margin_over_2se": bool(row[best] - row[run] > 2 * se) if se == se else False,
                     "same_winner_earlier_checkpoint": near.get(-a.checkpoint_step) == best, "same_winner_later_checkpoint": near.get(a.checkpoint_step) == best,
                     **{f"mcc_{l}": row[l] for l in row.index}})
    return pd.DataFrame(rows)


def plot(cells, per_cell, e_by, set_name, out, a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    embs = sorted({e for e, _ in cells}, key=lambda x: per_cell[per_cell.embedding == x].epoch.mean())
    losses = [l for l in LOSSES if any(l == k[1] for k in cells)]
    fig, axes = plt.subplots(len(embs), len(losses), figsize=(3.3 * len(losses), 2.2 * len(embs)), sharex=True, squeeze=False)
    for i, emb in enumerate(embs):
        for j, loss in enumerate(losses):
            ax = axes[i][j]
            epochs, mat = cells[(emb, loss)]
            vals = mat.mean(axis=1)
            grid, y = smooth_curve(epochs, vals, a.step, a.smooth, int(epochs.max()))
            ax.plot(grid, y / y.min(), color="0.35", lw=1)
            r = per_cell[(per_cell.embedding == emb) & (per_cell.loss_type == loss)].iloc[0]
            ax.axvline(r.epoch, color="tab:green", lw=1.4)
            if abs(r.unrestricted_argmin - r.epoch) > 50:
                ax.axvline(r.unrestricted_argmin, color="tab:red", ls="--", lw=1)
            ax.axvline(e_by[emb], color="tab:blue", ls=":", lw=1.2)
            ax.axvspan(a.max_search, grid.max(), color="0.9", zorder=0)
            ax.set_xscale("log"); ax.set_ylim(0.95, 3.0)
            if i == 0:
                ax.set_title(loss, fontsize=9)
            if j == 0:
                ax.set_ylabel(emb, fontsize=8)
    fig.suptitle(f"{set_name}: validation loss / its minimum.  green = epoch picked (before overfitting), red dashed = minimum of the whole run (late dip, ignored), "
                 f"blue = the epoch of this embedding (the average of its four green lines, moved to an MCC checkpoint); grey = beyond --max-search", fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=str(HERE / "train_val_loss_output"))
    ap.add_argument("--set", default="sam_full", help="file prefix of the sweep: sam_full, up ...")
    ap.add_argument("--out", default=str(HERE / "epoch_loss_selection_output"))
    ap.add_argument("--tol", type=float, default=0.03, help="relative rise above the running minimum that counts as overfitting")
    ap.add_argument("--patience", type=int, default=100, help="epochs the rise must last")
    ap.add_argument("--max-search", type=int, default=3000, help="hard limit of the search window (epochs)")
    ap.add_argument("--smooth", type=int, default=7, help="moving-average width in grid points")
    ap.add_argument("--step", type=int, default=10, help="grid step of the smoothed curve (epochs)")
    ap.add_argument("--aggregate", choices=["mean", "median"], default="mean")
    ap.add_argument("--epoch-scope", choices=["embedding", "global"], default="embedding",
                    help="embedding: every embedding gets its own epoch (the average over its losses of the epoch where the validation loss bottoms out); global: one epoch for all")
    ap.add_argument("--round", choices=["up", "nearest"], default="up",
                    help="how the averaged epoch becomes an MCC checkpoint: up = the next multiple of --round-to (153 -> 200; the default, because the MCC peaks later than the validation-loss minimum), nearest = the nearest checkpoint")
    ap.add_argument("--round-to", type=int, default=100, help="the multiple --round up goes to (use a multiple of --checkpoint-step so it is a checkpoint)")
    ap.add_argument("--epoch", type=float, help="use this epoch for every embedding instead of the computed ones (moved to the nearest checkpoint)")
    ap.add_argument("--loss-scope", choices=["cell", "hierarchical", "label", "global"], default="cell")
    ap.add_argument("--min-gain-se", type=float, default=2.0, help="hierarchical: standard errors a cell's own best loss must beat the global loss by")
    ap.add_argument("--mcc-window", type=int, default=0, help="average the MCC over the checkpoints within +-this many epochs of the chosen epoch")
    ap.add_argument("--bootstrap", type=int, default=200)
    ap.add_argument("--checkpoint-step", type=int, default=50)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--splits", type=int, default=0, help="select on the first N splits only (0 = all). The bake-off then keeps splits N+1.. as an unbiased check: nothing about the epoch or the losses was chosen on them")
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args()

    d = Path(a.dir)
    cf = sorted(d.glob(f"{a.set}__*__curves.csv.gz"))
    kf = sorted(d.glob(f"{a.set}__*__ckpt_mcc.csv.gz"))
    if not cf or not kf:
        sys.exit(f"no {a.set}__*__curves.csv.gz / __ckpt_mcc.csv.gz in {d}: run train_val_loss_curves.py first")
    curves = pd.concat([pd.read_csv(f, usecols=["embedding", "loss_type", "iteration", "epoch", "val_loss"]) for f in cf], ignore_index=True)
    ck = pd.concat([pd.read_csv(f, usecols=["embedding", "loss_type", "iteration", "epoch", "label", "mcc"]) for f in kf], ignore_index=True)
    if a.splits:
        curves, ck = curves[curves.iteration <= a.splits], ck[ck.iteration <= a.splits]
    grid = sorted(ck.epoch.unique())
    n_it = curves.iteration.nunique()
    cells = build_cells(curves)
    say(f"{a.set}: {len(cells)} (embedding, loss) cells, {n_it} splits" + (f" (the first {a.splits} of the sweep; the rest are left out of the selection)" if a.splits else "") + f", MCC checkpoints every {int(np.median(np.diff([g for g in grid if g <= 1000])))} epochs up to epoch 1000 ({len(grid)} checkpoints, {int(grid[0])}..{int(grid[-1])})")

    # ---- the epoch of every cell
    rows = []
    for (emb, loss), (epochs, mat) in cells.items():
        ep, found = cell_epoch(epochs, mat, np.arange(mat.shape[1]), a)
        fg, fy = smooth_curve(epochs, mat.mean(axis=1), a.step, a.smooth, int(epochs.max()))
        rows.append({"embedding": emb, "loss_type": loss, "epoch": ep, "onset_found": found, "unrestricted_argmin": float(fg[int(np.argmin(fy))])})
    per_cell = pd.DataFrame(rows)
    e_raw_by, e_by = epochs_by_embedding(per_cell, a.aggregate, grid, a.epoch_scope, a.epoch, a.round, a.round_to)
    rnd = (lambda v: snap(v, grid)) if a.epoch is not None else (lambda v: round_epoch(v, a.round, a.round_to, grid))
    e_raw = aggregate(list(e_raw_by.values()), a.aggregate)
    e_star = rnd(e_raw)                                         # one summary number for the whole sweep (the epoch itself when all embeddings share it)
    how_round = (f"rounded UP to the next multiple of {a.round_to}" if a.round == "up" and a.round_to else "moved to the nearest MCC checkpoint") if a.epoch is None else "moved to the nearest MCC checkpoint"
    say(f"\n1. EPOCH: overfitting onset = more than {a.tol:.0%} above the running minimum for {a.patience} epochs; search limited to the first {a.max_search} epochs")
    say(f"   epoch of the minimum before the onset, per cell: min {per_cell.epoch.min():.0f}, median {per_cell.epoch.median():.0f}, mean {per_cell.epoch.mean():.0f}, max {per_cell.epoch.max():.0f}")
    say("   by loss: " + ", ".join(f"{l} {per_cell[per_cell.loss_type == l].epoch.mean():.0f}" for l in LOSSES if (per_cell.loss_type == l).any()))
    say(f"   no overfitting onset found within the window in {int((~per_cell.onset_found).sum())} of {len(per_cell)} cells")
    late = per_cell[(per_cell.unrestricted_argmin - per_cell.epoch).abs() > 50]
    say(f"   the minimum of the WHOLE run differs from the pre-onset one in {len(late)} of {len(per_cell)} cells (a late dip or still-falling loss; ignored by the rule)"
        + (": " + ", ".join(f"{r.embedding}/{r.loss_type} {r.epoch:.0f} vs {r.unrestricted_argmin:.0f}" for r in late.head(6).itertuples()) if len(late) else ""))
    rng = np.random.default_rng(a.seed)
    boots = {e: [] for e in e_by}
    for _ in range(a.bootstrap):
        cols = rng.integers(0, n_it, n_it)
        cell_eps = {k: cell_epoch(ep, m, cols, a)[0] for k, (ep, m) in cells.items()}
        if a.epoch is not None:
            draw = dict(e_by)
        elif a.epoch_scope == "embedding":
            draw = {e: rnd(aggregate([v for (em, _), v in cell_eps.items() if em == e], a.aggregate)) for e in e_by}
        else:
            g0 = rnd(aggregate(list(cell_eps.values()), a.aggregate))
            draw = {e: g0 for e in e_by}
        for e, v in draw.items():
            boots[e].append(v)
    boots = {e: np.array(v) for e, v in boots.items()}
    if a.epoch is not None or a.epoch_scope == "global":
        e1 = next(iter(e_by.values()))
        say(f"   epoch for every embedding = {a.aggregate} of the cell epochs = {e_raw:.1f}" + (f"; set by hand to {a.epoch}" if a.epoch is not None else "") + f" -> {how_round}: checkpoint {e1:.0f}")
        b1 = next(iter(boots.values()))
        say(f"   bootstrap over the {n_it} splits ({a.bootstrap} draws): checkpoint {np.median(b1):.0f}, 95% interval {np.percentile(b1, 2.5):.0f}-{np.percentile(b1, 97.5):.0f}; "
            f"{(b1 == e1).mean():.0%} of the draws give {e1:.0f}")
    else:
        say(f"   the epoch of each EMBEDDING = {a.aggregate} of the epochs of its own losses, {how_round} (bootstrap over the {n_it} splits, {a.bootstrap} draws):")
        for e in sorted(e_by, key=lambda x: e_by[x]):
            b = boots[e]
            say(f"      {e:16s} {e_raw_by[e]:6.1f} -> checkpoint {e_by[e]:.0f} | 95% interval {np.percentile(b, 2.5):.0f}-{np.percentile(b, 97.5):.0f} | {(b == e_by[e]).mean():.0%} of the draws give it")

    # ---- the loss
    sc = scores_per_embedding(ck, e_by, a.mcc_window)
    mapping = choose(sc, a.loss_scope, a.min_gain_se)
    table = cell_table(sc, mapping, ck, e_by, a.mcc_window, grid, a)
    sel_in = score_map(sc, mapping)
    held, fixed, avg = leave_one_split_out(sc, a.loss_scope, a.min_gain_se)
    say("\n2. LOSS at " + (f"epoch {e_star:.0f}" if len(set(e_by.values())) == 1 else "each embedding's own epoch") + (f" (MCC averaged over +-{a.mcc_window} epochs)" if a.mcc_window else "") + f", scope '{a.loss_scope}'" + (f", min gain {a.min_gain_se} SE" if a.loss_scope == 'hierarchical' else ""))
    say("   chosen losses: " + ", ".join(f"{l} {sum(1 for v in mapping.values() if v == l)}" for l in LOSSES))
    g = choose(sc, "global", a.min_gain_se)
    say(f"   best single loss for everything at these epochs: {next(iter(g.values()))} (mean MCC {score_map(sc, g):.4f})")
    say(f"   the winner's margin over the runner-up: median {table.margin.median():.4f}; paired standard error over the splits: median {table.paired_se.median():.4f}; "
        f"margin above 2 standard errors in {int(table.margin_over_2se.sum())} of {len(table)} cells")
    say(f"   the winner is the same one checkpoint earlier in {int(table.same_winner_earlier_checkpoint.sum())} cells, one checkpoint later in {int(table.same_winner_later_checkpoint.sum())} (of {len(table)})")
    say(f"   macro MCC of the chosen map: in-sample {sel_in:.4f}")
    say(f"   leave-one-split-out (choose on the other splits, score on the held-out one): chosen map {held:.4f} | best single loss {fixed:.4f} | average over the losses {avg:.4f}")
    say(f"   -> optimism of the in-sample score: {sel_in - held:+.4f}. Report the held-out number, or select on one set of splits and report on another.")
    if a.loss_scope in ("cell", "label") and fixed - held > 0.002:  # a smaller gap is noise, not a reason to warn
        say(f"   NOTE: one loss for everything scores higher on held-out splits ({fixed:.4f}) than this per-{a.loss_scope} map ({held:.4f}): choosing a loss per {a.loss_scope} is fitting noise "
            f"at this number of splits. --loss-scope hierarchical keeps the single best loss and changes it only where a cell's gain is real.")

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    loss_map = {}
    for (emb, lab), loss in mapping.items():
        loss_map.setdefault(emb, {})[lab] = loss
    (out / f"loss_map_{a.set}.json").write_text(json.dumps(loss_map, indent=1))
    (out / f"selection_{a.set}.json").write_text(json.dumps({
        "set": a.set, "epoch": e_star, "epoch_unsnapped": e_raw, "epoch_scope": "global" if (a.epoch is not None or a.epoch_scope == "global") else "embedding",
        "epoch_by_embedding": {e: float(v) for e, v in e_by.items()}, "epoch_unsnapped_by_embedding": {e: float(v) for e, v in e_raw_by.items()}, "aggregate": a.aggregate,
        "epoch_rounding": "manual (nearest checkpoint)" if a.epoch is not None else a.round, "epoch_round_to": a.round_to if a.epoch is None and a.round == "up" else None,
        "epoch_interval_95_by_embedding": {e: [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))] for e, b in boots.items()},
        "selected_on_first_n_splits": a.splits or None, "n_splits_used": int(n_it), "loss_scope": a.loss_scope, "min_gain_se": a.min_gain_se, "mcc_window": a.mcc_window, "overfitting_tolerance": a.tol, "patience": a.patience, "max_search": a.max_search,
        "in_sample_macro_mcc": sel_in, "held_out_macro_mcc": held, "loss_map": loss_map}, indent=1))
    table.to_csv(out / f"selection_cells_{a.set}.csv", index=False)
    per_cell.to_csv(out / f"epoch_per_cell_{a.set}.csv", index=False)
    (out / f"selection_report_{a.set}.txt").write_text("\n".join(LINES) + "\n")
    if not a.no_plot:
        plot(cells, per_cell, e_by, a.set, out / f"epoch_selection_{a.set}.png", a)
    say(f"\nwrote {out}/ (selection_{a.set}.json, loss_map_{a.set}.json, selection_cells_{a.set}.csv, epoch_per_cell_{a.set}.csv, report, figure)")


if __name__ == "__main__":
    main()
