#!/usr/bin/env python3
"""The dark proteins of a release, as the inputs of CLEAN and Horizyn.

  python make_dark_inputs.py --release <release> [--name dark53]

Reads   <release>/umap/dark_clusters.csv                    the dark clusters (5_umap/cluster_selection.py): cluster, lineage
        <release>/umap/sam_mtase_umap_prott5_coords.csv     every protein of the map with its cluster (Entry, Organism, cluster)
        <release>/embedding_inputs/proteomes_ipr029063.fasta   the sequences that were embedded
Writes  <release>/dark_annotation/<name>_members.csv        Entry, Organism, cluster, lineage   (rows ordered by dark_clusters.csv, then as on the map)
        <release>/dark_annotation/<name>.fasta              one single-line record per member
<name> defaults to dark<number of dark clusters>, e.g. dark53: the file stem CLEAN uses and the prefix of its result files (run_tables.sh reads the same name).
"""
import argparse
from pathlib import Path

import pandas as pd


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for ln in open(path):
        ln = ln.rstrip("\n")
        if ln.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name, buf = ln[1:].split()[0], []
        else:
            buf.append(ln.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def build(dark, coords, seqs):
    """dark: dark_clusters.csv; coords: the map's table; seqs: {Entry: sequence} -> members DataFrame, sequences of the members."""
    order = {int(c): i for i, c in enumerate(dark["cluster"])}
    lineage = dict(zip(dark["cluster"].astype(int), dark["lineage"]))
    m = coords[coords["cluster"].isin(order)][["Entry", "Organism", "cluster"]].copy()
    m["cluster"] = m["cluster"].astype(int)
    m["lineage"] = m["cluster"].map(lineage)
    m["_o"] = m["cluster"].map(order)
    m = m.sort_values("_o", kind="stable").drop(columns="_o").reset_index(drop=True)
    missing = [e for e in m["Entry"] if e not in seqs]
    if missing:
        raise SystemExit(f"{len(missing)} dark proteins have no sequence in the embedding inputs, e.g. {missing[:3]}")
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", required=True)
    ap.add_argument("--name", default=None)
    a = ap.parse_args()
    rel = Path(a.release)
    dark = pd.read_csv(rel / "umap" / "dark_clusters.csv")
    if dark.empty:
        raise SystemExit("no dark clusters in this release's map: nothing to annotate")
    coords = pd.read_csv(rel / "umap" / "sam_mtase_umap_prott5_coords.csv", usecols=["Entry", "Organism", "cluster"], low_memory=False)
    seqs = read_fasta(rel / "embedding_inputs" / "proteomes_ipr029063.fasta")
    m = build(dark, coords, seqs)
    name = a.name or f"dark{len(dark)}"
    out = rel / "dark_annotation"
    out.mkdir(exist_ok=True)
    m.to_csv(out / f"{name}_members.csv", index=False)
    with open(out / f"{name}.fasta", "w") as f:
        for e in m["Entry"]:
            f.write(f">{e}\n{seqs[e]}\n")
    (out / "dark_name.txt").write_text(name + "\n")      # the stem the later steps (CLEAN, Horizyn, tables) read
    print(f"{name}: {len(dark)} dark clusters ({dark['lineage'].value_counts().to_dict()}), {len(m):,} proteins -> {out}")


if __name__ == "__main__":
    main()
