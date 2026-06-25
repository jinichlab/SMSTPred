#!/usr/bin/env python3
"""
Global ProstT5 UMAP of all tree-of-life SAM-MTases, colored by substrate type.

ProstT5 counterpart to sam_mtase_umap_2026-05-26 (which used ProtT5). Loads every
mean-pooled ProstT5 vector from the tree-of-life embedding set
(reference_files/prostt5_embeddings/{organism}_embeddings.h5), assigns each
protein a substrate from the Stage-3 ProstT5 deployment calls (multi-label ->
single label / "Mixed" / "Unassigned"), runs the SAME UMAP+HDBSCAN params as the
ProstT5 central re-embedding, and writes a static PNG
HTML, and a coords CSV. Both embedding AND labels are ProstT5, so the figure is
internally consistent with the new classifier.

Outputs:
  sam_mtase_umap_prostt5_coords.csv
  sam_mtase_umap_prostt5_static.png
  sam_mtase_umap_prostt5_interactive.html

Env: UMAP (umap-learn). Interactive HTML: substrate_umap_datamap.py.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
import umap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.cluster import HDBSCAN

EMB_DIR = Path("/home/nroethler/Code/classifier/reference_files/prostt5_embeddings")
PRED_DIR = Path("/home/nroethler/Code/journal_article/repo/3_predictions/predictions_output")
META_CSV = Path("/home/nroethler/Code/classifier/reference_files/selected_1259_organisms.csv")
IO = Path(__file__).parent.resolve()
EMB_SUFFIX = "_embeddings.h5"
PRED_SUFFIX = "_embeddings_predictions.csv"

RANDOM_SEED, N_NEIGHBORS, MIN_DIST, METRIC = 42, 15, 0.1, "cosine"
HDBSCAN_MIN_CLUSTER_SIZE, HDBSCAN_MIN_SAMPLES = 50, 10
LABELS = ["Other", "Protein", "RNA", "DNA"]
SUBSTRATE_COLORS = {
    "RNA": "#cb181d", "DNA": "#6a51a3", "Protein": "#238b45",
    "Other": "#e6550d", "Mixed": "#636363", "Unassigned": "#dddddd",
}


def organism_to_stem(name):
    return str(name).replace(" ", "_").replace("/", "_")


def derive_substrate(row):
    pos = [l for l in LABELS if row[f"call_{l}"] == 1]
    return "Unassigned" if not pos else (pos[0] if len(pos) == 1 else "Mixed")


def load_predictions():
    meta = pd.read_csv(META_CSV)
    stem_to_org = {organism_to_stem(o): o for o in meta["Organism"]}
    frames = []
    for f in sorted(PRED_DIR.glob(f"*{PRED_SUFFIX}")):
        d = pd.read_csv(f, low_memory=False)
        if d.empty or not {f"call_{l}" for l in LABELS}.issubset(d.columns):
            continue
        d["Organism"] = stem_to_org.get(f.name[:-len(PRED_SUFFIX)], "")
        d["Substrate"] = d.apply(derive_substrate, axis=1)
        frames.append(d[["Entry", "Organism", "Substrate",
                         "score_RNA", "score_DNA", "score_Protein", "score_Other"]])
    pred = pd.concat(frames, ignore_index=True).drop_duplicates("Entry").set_index("Entry")
    print(f"predictions: {len(pred):,} proteins | "
          f"{pred['Substrate'].value_counts().to_dict()}", flush=True)
    return pred


def load_embeddings(keep_entries):
    ents, vecs = [], []
    for h in sorted(EMB_DIR.glob(f"*{EMB_SUFFIX}")):
        with h5py.File(h, "r") as f:
            for sid in f:
                if sid not in keep_entries:
                    continue
                g = f[sid]
                v = np.array(g["mean"]) if "mean" in g else np.array(g)
                if v.ndim == 1 and v.shape[0] == 1024 and np.isfinite(v).all():
                    ents.append(sid); vecs.append(v.astype(np.float32))
    print(f"embeddings: {len(ents):,} loaded", flush=True)
    return ents, np.vstack(vecs)


def main():
    pred = load_predictions()
    ents, X = load_embeddings(set(pred.index))
    df = pd.DataFrame({"Entry": ents})
    df["Organism"] = pred.loc[df.Entry, "Organism"].values
    df["Substrate"] = pred.loc[df.Entry, "Substrate"].values
    for l in LABELS:
        df[f"score_{l}"] = pred.loc[df.Entry, f"score_{l}"].values
    print(f"aligned matrix: {X.shape} -> UMAP", flush=True)

    coords = umap.UMAP(n_components=2, n_neighbors=N_NEIGHBORS, min_dist=MIN_DIST,
                       metric=METRIC, random_state=RANDOM_SEED, verbose=True).fit_transform(X)
    df["UMAP_1"], df["UMAP_2"] = coords[:, 0], coords[:, 1]
    df["cluster"] = HDBSCAN(min_cluster_size=HDBSCAN_MIN_CLUSTER_SIZE,
                            min_samples=HDBSCAN_MIN_SAMPLES).fit_predict(coords)
    nclu = df.cluster.nunique() - (1 if (df.cluster == -1).any() else 0)
    print(f"{nclu} clusters, {int((df.cluster == -1).sum()):,} noise", flush=True)
    df.to_csv(IO / "sam_mtase_umap_prostt5_coords.csv", index=False)

    # --- static PNG ---
    fig, ax = plt.subplots(figsize=(14, 11))
    order = ["Unassigned", "Mixed", "Other", "Protein", "DNA", "RNA"]
    dim = {"Unassigned": (0.20, 4), "Mixed": (0.40, 6)}
    for k, sub in enumerate(order):
        m = df["Substrate"] == sub
        if not m.any():
            continue
        a, s = dim.get(sub, (0.7, 7))
        ax.scatter(df.loc[m, "UMAP_1"], df.loc[m, "UMAP_2"], c=SUBSTRATE_COLORS[sub],
                   label=f"{sub} (n={int(m.sum()):,})", alpha=a, s=s,
                   edgecolors="none", zorder=2 + k)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.set_title(f"Tree-of-life SAM-MTases — ProstT5 (structure-aware)\n"
                 f"{len(df):,} proteins, colored by ProstT5 substrate call",
                 fontweight="bold")
    ax.legend(loc="upper right", fontsize=10, framealpha=0.9)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    plt.tight_layout()
    plt.savefig(IO / "sam_mtase_umap_prostt5_static.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("wrote static PNG + coords CSV", flush=True)
    # The interactive figure is rendered separately by substrate_umap_datamap.py
    # (datamapplot), which reads this coords CSV.


if __name__ == "__main__":
    main()
