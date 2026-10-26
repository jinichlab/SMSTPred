#!/usr/bin/env python3
"""
Per-label MCC vs epoch grids for the MLP substrate classifiers (one figure family).

Writes two consistently-styled grids (plus a third with --balanced):
  mlp_embedding_comparison_output/full_mcc_curves_grid.png       SAM full
  mlp_embedding_comparison_output/balanced_mcc_curves_grid.png   SAM balanced (only with --balanced; not in the manuscript)
        -> rows = embedding (ProtT5 / ProstT5 / ESM-IF / the two concatenations),
           cols = label, lines = the 4 loss types
  UP_prostt5_loss_comparison_output/UP_prostt5_mcc_curves_grid.png  UP (ProstT5)
        -> single ProstT5 row, cols = label, lines = loss

Each panel: MCC vs epoch (log x), mean over the cluster-disjoint iterations with a
+/-1 SD band, and a dot at every loss's peak MCC. Pure CSV -> PNG, no GPU.

Sources:
  mlp_embedding_comparison_output/mlp_embedding_comparison_per_label.csv
      (SAM embedding sweep; cols: dataset, embedding, loss_type, iteration, epochs,
       label, mcc, ...)
  UP_prostt5_loss_comparison_output/UP_prostt5_mlp_loss_comparison_per_label.csv
      (UP ProstT5 loss sweep; single embedding -> embedding column added here)
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent
LABELS = ["DNA", "RNA", "Protein", "Other"]
LOSSES = ["bce", "weighted_bce", "sqrt_weighted_bce", "focal"]
LOSS_COLORS = dict(zip(LOSSES, ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]))


def render_grid(pl, embeddings, outpath, ylim=None):
    n_iters = pl["iteration"].nunique()
    nr, nc = len(embeddings), len(LABELS)
    fig, axes = plt.subplots(nr, nc, figsize=(4.4 * nc, 3.0 * nr + 1.0),
                             squeeze=False, sharex=True, sharey=(ylim is not None))
    for r, emb in enumerate(embeddings):
        for c, lab in enumerate(LABELS):
            ax = axes[r][c]
            sub = pl[(pl.embedding == emb) & (pl.label == lab)]
            for loss in LOSSES:
                g = sub[sub.loss_type == loss].groupby("epochs").mcc.agg(["mean", "std"])
                if g.empty:
                    continue
                x = g.index.values
                m = g["mean"].values
                s = g["std"].fillna(0).values
                col = LOSS_COLORS[loss]
                ax.plot(x, m, color=col, lw=1.7, label=loss, zorder=3)
                ax.fill_between(x, m - s, m + s, color=col, alpha=0.15, lw=0, zorder=1)
                pk = m.argmax()
                ax.plot(x[pk], m[pk], "o", color=col, ms=5, zorder=4)  # peak
            ax.set_xscale("log")
            if ylim is not None:
                ax.set_ylim(*ylim)
            ax.grid(alpha=0.3, which="both")
            if r == 0:
                ax.set_title(lab, fontsize=12, fontweight="bold")
            if c == 0:
                ax.set_ylabel(f"{emb}\nMCC", fontsize=9)
            if r == nr - 1:
                ax.set_xlabel("epoch (log)")
    axes[0][-1].legend(fontsize=7, loc="lower right", title=f"loss (n={n_iters})")
    fig.tight_layout()
    fig.savefig(outpath, dpi=140)
    plt.close(fig)
    print("wrote", outpath)


def main():
    # --- SAM embedding sweep: rows = embedding, sharey to compare embeddings ---
    emb_dir = HERE / "mlp_embedding_comparison_output"
    epl = pd.read_csv(emb_dir / "mlp_embedding_comparison_per_label.csv")
    # SaProt was swept separately (mlp_embedding_comparison_saprot.py -> mlp_embedding_comparison_saprot_output/);
    # identical schema/epochs/iterations -> fold it in as one more embedding row.
    sap_csv = HERE / "mlp_embedding_comparison_saprot_output" / "mlp_embedding_comparison_per_label.csv"
    if sap_csv.exists():
        epl = pd.concat([epl, pd.read_csv(sap_csv)[epl.columns]], ignore_index=True)
    EMBS = ["ProstT5", "SaProt", "ProtT5", "ESM-IF", "ProtT5+ESM-IF", "ProstT5+ESM-IF"]
    for ds in (["full", "balanced"] if "--balanced" in sys.argv else ["full"]):
        render_grid(epl[epl.dataset == ds], EMBS,
                    emb_dir / f"{ds}_mcc_curves_grid.png", ylim=(0.0, 1.0))

    # --- UP ProstT5 loss sweep: single embedding, per-panel y to show loss split ---
    up_dir = HERE / "UP_prostt5_loss_comparison_output"
    upl = pd.read_csv(up_dir / "UP_prostt5_mlp_loss_comparison_per_label.csv")
    upl["embedding"] = "ProstT5"
    render_grid(upl, ["ProstT5"],
                up_dir / "UP_prostt5_mcc_curves_grid.png", ylim=None)


if __name__ == "__main__":
    main()
