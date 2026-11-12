#!/usr/bin/env python3
"""Global UMAP of every tree-of-life SAM-MTase of a refresh release, from its ProtT5 embedding, coloured by the ensemble's substrate call.

  python global_umap.py [--release <folder>] [--embedding ProtT5] [--frozen-map <coords csv>] [--recompute]      (default release: $RELEASE_DIR, else the newest release_* of $UNIPROT_PULL_DIR)

THE MAP IS FROZEN BY DEFAULT. UMAP + HDBSCAN are not reproducible from run to run (a rerun on the same inputs, seed and environment moves the points by tens of units and gives 933-947
clusters), and everything after them (dark clusters, the Fig 8 clusters, motif counts) depends on the exact map. So when <release>/umap/sam_mtase_umap_<emb>_coords.csv (or --frozen-map) exists
and covers exactly the proteins of this release, its UMAP_1 / UMAP_2 / cluster columns are reused as they are (only the scores and calls are refreshed from the current predictions) and no UMAP
is run. A release whose proteins differ from the map's, or no map at all, gets a new UMAP (a new release gets a new map); --recompute (or UMAP_RECOMPUTE=1) forces a new one.

The release counterpart of global_reumap_prostt5.py, with the SAME settings (UMAP: 15 neighbours, min_dist 0.1, cosine, seed 42; HDBSCAN on the 2-D coordinates:
min_cluster_size 50, min_samples 10). Reads <release>/embeddings/embeddings_proteomes_<embedding>.h5 (mean-pooled 1,024-d vector per Entry, keyed by Entry) and
<release>/predictions_all_proteomes.csv.gz (the ProtT5 ensemble: ALL proteins <= 2,000 aa with their calls, in domain or not). Substrate = the single class called;
"Mixed" = two or more; "None" = no call. Writes to <release>/umap/:
  sam_mtase_umap_<emb>_coords.csv           Entry, Organism, proteome, Substrate, the four scores and calls, disagreement, in_domain, UMAP_1/2, HDBSCAN cluster
  sam_mtase_umap_<emb>_static_datamap.png   the static map (same style as the archived figure)
  sam_mtase_umap_<emb>_interactive.html     the datamapplot map (hover, search, click-to-UniProt, legend toggles)
  umap_report.txt                           counts, clusters, and how well the 2-D map keeps the substrate classes together (neighbour purity)
"""
import json
import argparse
import os
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import umap
from sklearn.cluster import HDBSCAN
from sklearn.neighbors import NearestNeighbors

import substrate_umap_datamap as sd

REPO = Path(__file__).resolve().parents[1]
META_CSV = Path(os.environ.get("SMSTPRED_DATA", REPO / "data")) / "selected_1259_organisms.csv"
RANDOM_SEED, N_NEIGHBORS, MIN_DIST, METRIC = 42, 15, 0.1, "cosine"
HDBSCAN_MIN_CLUSTER_SIZE, HDBSCAN_MIN_SAMPLES = 50, 10
LABELS = ["Metabolite", "Protein", "RNA", "DNA"]


def substrate_of(calls):
    """calls: (n, 4) 0/1 array in LABELS order -> Substrate strings."""
    n = calls.sum(axis=1)
    single = np.array(LABELS)[calls.argmax(axis=1)]
    return np.where(n == 0, "None", np.where(n == 1, single, "Mixed"))


def load(rel, embedding):
    pred = pd.read_csv(rel / "predictions_all_proteomes.csv.gz")
    summ = pd.read_csv(rel / "proteome_summary.csv")[["proteome_id", "stem"]]
    meta = pd.read_csv(META_CSV)[["Proteome Id", "Organism"]].merge(summ, left_on="Proteome Id", right_on="proteome_id")
    pred = pred.merge(meta[["stem", "Organism"]], left_on="proteome", right_on="stem", how="left").drop(columns="stem")
    pred["Organism"] = pred["Organism"].fillna(pred["proteome"])
    h5 = rel / "embeddings" / f"embeddings_proteomes_{embedding}.h5"
    with h5py.File(h5, "r") as f:
        have = [e for e in pred.Entry if e in f]
        X = np.vstack([np.asarray(f[e], dtype=np.float32) for e in have])
    pred = pred.set_index("Entry").loc[have].reset_index()
    ok = np.isfinite(X).all(axis=1)
    print(f"{len(pred):,} predictions, {len(have):,} with an embedding, {int((~ok).sum())} non-finite dropped", flush=True)
    return pred[ok].reset_index(drop=True), X[ok]


def release_label(rel):
    """release_<UniProt release> from the release folder's pull_manifest.json (the folder itself may be a scratch copy with another name)."""
    try:
        return "release_" + json.loads((rel / "pull_manifest.json").read_text())["uniprot_release"]
    except Exception:
        return rel.name


def neighbor_purity(coords, substrate, k=15):
    """Of the points with a single-class call, the share whose k nearest neighbours in the 2-D map carry the same class (mean over points), and the same for each class."""
    m = np.isin(substrate, LABELS)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords[m])
    idx = nn.kneighbors(coords[m], return_distance=False)[:, 1:]
    lab = substrate[m]
    same = (lab[idx] == lab[:, None]).mean(axis=1)
    return float(same.mean()), {l: float(same[lab == l].mean()) for l in LABELS if (lab == l).any()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release")
    ap.add_argument("--embedding", default="ProtT5")
    ap.add_argument("--frozen-map", help="coordinates and clusters to reuse (default: <release>/umap/sam_mtase_umap_<embedding>_coords.csv when it exists)")
    ap.add_argument("--recompute", action="store_true", help="compute a new UMAP and HDBSCAN clustering even when a map exists (also UMAP_RECOMPUTE=1)")
    a = ap.parse_args()
    sys.path.insert(0, str(REPO / "1_clustering"))
    from mmseqs_clusters import newest_release
    rel = Path(a.release or os.environ.get("RELEASE_DIR") or newest_release(Path(os.environ.get("UNIPROT_PULL_DIR", REPO / "data" / "uniprot_pull")), os.environ.get("RELEASE")))
    out = rel / "umap"
    out.mkdir(exist_ok=True)
    emb = a.embedding.lower()

    df, X = load(rel, a.embedding)
    calls = df[[f"call_{l}" for l in LABELS]].to_numpy().astype(int)
    df["Substrate"] = substrate_of(calls)
    print("substrate:", df.Substrate.value_counts().to_dict(), flush=True)

    frozen = None
    fpath = Path(a.frozen_map) if a.frozen_map else out / f"sam_mtase_umap_{emb}_coords.csv"
    if not (a.recompute or os.environ.get("UMAP_RECOMPUTE") == "1") and fpath.exists():
        fm = pd.read_csv(fpath, usecols=["Entry", "UMAP_1", "UMAP_2", "cluster"]).drop_duplicates("Entry").set_index("Entry")
        if set(fm.index) == set(df.Entry):
            frozen = fm.loc[df.Entry]
            print(f"FROZEN MAP: coordinates and clusters of {len(fm):,} proteins reused from {fpath} (no UMAP run; --recompute makes a new map)", flush=True)
        else:
            print(f"the map in {fpath} covers {len(fm):,} proteins, this release has {len(df):,} (only {len(set(fm.index) & set(df.Entry)):,} in common): computing a new map", flush=True)
    t0 = time.time()
    if frozen is not None:
        coords = frozen[["UMAP_1", "UMAP_2"]].to_numpy(dtype=np.float32)
        df["UMAP_1"], df["UMAP_2"], df["cluster"] = coords[:, 0], coords[:, 1], frozen["cluster"].to_numpy()
        map_note = f"map: frozen, reused from {fpath.name} (not recomputed)"
    else:
        coords = umap.UMAP(n_components=2, n_neighbors=N_NEIGHBORS, min_dist=MIN_DIST, metric=METRIC, random_state=RANDOM_SEED, verbose=True).fit_transform(X)
        print(f"UMAP {time.time() - t0:.0f} s", flush=True)
        df["UMAP_1"], df["UMAP_2"] = coords[:, 0], coords[:, 1]
        df["cluster"] = HDBSCAN(min_cluster_size=HDBSCAN_MIN_CLUSTER_SIZE, min_samples=HDBSCAN_MIN_SAMPLES).fit_predict(coords)
        map_note = "map: computed in this run (UMAP and HDBSCAN are not reproducible from run to run: keep this map, or freeze it with the release)"
    nclu = df.cluster.nunique() - (1 if (df.cluster == -1).any() else 0)
    keep = ["Entry", "Organism", "proteome", "Substrate", *[f"score_{l}" for l in LABELS], *[f"call_{l}" for l in LABELS], "disagreement", "in_domain", "Length",
            "UMAP_1", "UMAP_2", "cluster"]
    df[keep].to_csv(out / f"sam_mtase_umap_{emb}_coords.csv", index=False)

    pur, per = neighbor_purity(coords, df.Substrate.to_numpy())
    lines = [f"{a.embedding} UMAP of {len(df):,} proteins from {df.proteome.nunique()} proteomes ({release_label(rel)})",
             f"UMAP: n_neighbors {N_NEIGHBORS}, min_dist {MIN_DIST}, metric {METRIC}, seed {RANDOM_SEED}; HDBSCAN min_cluster_size {HDBSCAN_MIN_CLUSTER_SIZE}, min_samples {HDBSCAN_MIN_SAMPLES}",
             map_note,
             f"substrate: {df.Substrate.value_counts().to_dict()}",
             f"HDBSCAN: {nclu} clusters, {int((df.cluster == -1).sum()):,} unclustered",
             f"neighbour purity (single-call proteins; share of the 15 nearest neighbours in the 2-D map with the same call): {pur:.3f}; " + ", ".join(f"{l} {v:.3f}" for l, v in per.items())]
    (out / "umap_report.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)

    # sd works on a "Membership" list per point (the deployment's calls); here the calls are in the table
    df["Membership"] = [[l for l, c in zip(LABELS, row) if c] or ["None"] for row in calls]
    sd.create_static_plot(df, out / f"sam_mtase_umap_{emb}_static_datamap.png")
    sd.create_interactive_plot(df, out / f"sam_mtase_umap_{emb}_interactive.html",
                               title=f"SAM Methyltransferases by Predicted Substrate Type ({a.embedding})",
                               sub_title=f"UMAP projection of {len(df):,} tree-of-life SAM MTases ({a.embedding} sequence embeddings)")
    print("done", flush=True)


if __name__ == "__main__":
    main()
