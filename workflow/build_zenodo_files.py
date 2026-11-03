#!/usr/bin/env python3
"""Assemble the files of the Zenodo data deposit from a release folder.

  python workflow/build_zenodo_files.py --release <release> --out <folder>      (default release: $RELEASE_DIR, else the newest; the folder must be empty or not exist)

Writes into <folder> exactly what is to be uploaded: the large inputs and results that do not belong in git (embeddings, labels, predictions, the deployed ensemble, the Horizyn scores of the dark
proteins, the dark-cluster tables, motif and structure results, benchmark tables, the map, the tree, the Rhea / ENZYME tables), a README.md that describes each file, and MANIFEST.sha256.
Symbolic links are never followed (the Horizyn run folder links into an external checkout), and nothing outside the release folder is read. The upload itself is zenodo_upload/zenodo_draft.py.
"""
import argparse
import hashlib
import json
import os
import shutil
import sys
import tarfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "1_clustering"))

# (output name, kind, source relative to the release | list of sources, description)
FILES = [
    ("embeddings_SAM_ProtT5.h5", "copy", "embeddings/embeddings_SAM_ProtT5.h5", "ProtT5 embeddings of the SAM-methyltransferase set (the deployed model's features), one 1,024-d vector per UniProt accession"),
    ("embeddings_SAM_ProstT5.h5", "copy", "embeddings/embeddings_SAM_ProstT5.h5", "ProstT5 embeddings of the same set"),
    ("embeddings_SAM_SaProt.h5", "copy", "embeddings/embeddings_SAM_SaProt.h5", "SaProt embeddings (structure-aware; from the structures in sam_structures_RELEASE.tar.gz)"),
    ("embeddings_SAM_ESM-IF.h5", "copy", "embeddings/embeddings_SAM_ESM-IF.h5", "ESM-IF embeddings (from the same structures)"),
    ("embeddings_UP_ProtT5.h5", "copy", "embeddings/embeddings_UP_ProtT5.h5", "ProtT5 embeddings of the broad UniProt-enzyme (UP) benchmark set"),
    ("embeddings_proteomes_ProtT5.h5", "copy", "embeddings/embeddings_proteomes_ProtT5.h5", "ProtT5 embeddings of the methyltransferase-domain proteins of the reference proteomes (keyed by accession; proteome_members.csv gives the proteome of each)"),
    ("proteome_members.csv", "copy", "embedding_inputs/proteome_members.csv", "proteome, accession and length of every protein in the proteome embeddings"),
    ("labels_SAM_RELEASE.csv", "copy", "labels_SAM.csv", "SAM-methyltransferase training set with substrate-type labels (Metabolite / Protein / RNA / DNA) typed from Rhea"),
    ("labels_UP_RELEASE.csv", "copy", "labels_UP.csv", "UP comparison set with its substrate-type labels"),
    ("predictions_all_proteomes.csv.gz", "copy", "predictions_all_proteomes.csv.gz", "the deployed ensemble's scores, calls, ensemble disagreement and applicability-domain flag for every protein of the proteome embeddings"),
    ("deployment_ensemble_RELEASE.tar", "tar", "deployment", "the deployed ProtT5 6-member ensemble: weights, calibration and applicability-domain reference"),
    ("horizyn_dark_scores.npz", "copy", "dark_annotation/horizyn_run/horizyn_full_scores.npz", "Horizyn scores of every dark-cluster protein against all Rhea reaction entries (scores, accessions, clusters, rhea_ids)"),
    ("horizyn_dark_top50.csv", "copy", "dark_annotation/horizyn_run/horizyn_full_topk.csv", "the 50 best-scoring reactions of each dark-cluster protein"),
    ("dark_annotation_RELEASE.tar", "tar", ["dark_annotation"], "dark-cluster members and FASTA, CLEAN results, Horizyn family tables, the notes and Tables S1 / S2 (without the Horizyn run folder)"),
    ("motif_analysis_RELEASE.tar", "tar", ["motif_analysis"], "STREME motifs, FIMO status of every protein, the specificity matrix, Table 2 and Table S3 inputs, the structure search and the structure picks"),
    ("sam_structures_RELEASE.tar.gz", "targz", "structures", "the structure store: one <Entry>.pdb per benchmark protein (AlphaFold DB models, ColabFold predictions) and structures_status.csv"),
    ("benchmark_RELEASE.tar", "tar", ["benchmark", "sam_vs_up"], "the bake-off and the MTase-vs-UP comparison: per-split tables, statistics, loss curves"),
    ("umap_RELEASE.tar", "tar", ["umap"], "the ProtT5 map (coordinates, HDBSCAN clusters, dark / metazoan-protein / vertebrate-DNA cluster tables)"),
    ("tree_RELEASE.tar", "tar", ["tree"], "the taxonomy tree (Newick) and the per-organism annotation tables of Fig 3"),
    ("reference_tables_RELEASE.tar", "tar", ["reference_tables"], "ExPASy ENZYME, Rhea EC links, Rhea direction table and Rhea equation text used to name enzymes and reactions"),
    ("pull_manifest.json", "copy", "pull_manifest.json", "what was pulled from UniProt and when"),
]
SKIP_DIRS = {"horizyn_run", "cif", "__pycache__", "tmp_mmseqs"}
SKIP_SUFFIX = (".part",)
SKIP_NAMES = {"context_cache.json"}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def add_tree(tar, rel, src):
    """Add folder `src` (relative to the release) without symlinks, caches and the big regenerable files."""
    for root, dirs, files in os.walk(rel / src):
        dirs[:] = sorted(d for d in dirs if d not in SKIP_DIRS and not (Path(root) / d).is_symlink())
        for f in sorted(files):
            p = Path(root) / f
            if p.is_symlink() or f in SKIP_NAMES or f.endswith(SKIP_SUFFIX):
                continue
            if f in ("corpus.csv", "corpus_all.fasta", "corpus_every.clean.fasta") or f.startswith(("class_", "reps_", "ctrl_", "clust_")):
                continue         # regenerable from the pull and large
            tar.add(p, arcname=str(p.relative_to(rel)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    from mmseqs_clusters import newest_release
    rel = Path(a.release or os.environ.get("RELEASE_DIR") or newest_release(Path(os.environ.get("UNIPROT_PULL_DIR", REPO / "data" / "uniprot_pull")), os.environ.get("RELEASE")))
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit(f"{out} is not empty")
    out.mkdir(parents=True, exist_ok=True)
    tag = rel.name.replace("release_", "")
    man = json.load(open(rel / "pull_manifest.json"))
    rows = []
    for name, kind, src, desc in FILES:
        name, desc = name.replace("RELEASE", tag), desc.replace("RELEASE", tag)
        dst = out / name
        if kind == "copy":
            s = rel / src
            if not s.exists():
                print(f"  skip {name}: {s} not found")
                continue
            try:
                os.link(s, dst)                  # no copy when it is the same file system
            except OSError:
                shutil.copyfile(s, dst)
        elif kind in ("tar", "targz"):
            srcs = src if isinstance(src, list) else [src]
            if not all((rel / x).exists() for x in srcs):
                print(f"  skip {name}: {[x for x in srcs if not (rel / x).exists()]} not found")
                continue
            with tarfile.open(dst, "w:gz" if kind == "targz" else "w") as t:
                for x in srcs:
                    add_tree(t, rel, x)
        rows.append((name, dst.stat().st_size, desc))
        print(f"  {name}  {dst.stat().st_size / 1e6:9.1f} MB", flush=True)
    readme = [f"# SMSTPred data, UniProt release {man['uniprot_release']} ({man.get('uniprot_release_date', '')})", "",
              "Data accompanying the manuscript on SAM-methyltransferase substrate-type prediction (N. Roethler, A. Jinich): precomputed protein language-model embeddings, the substrate-type label tables,",
              "the deployed model, its predictions over the reference proteomes, and the results of the analyses built on them. Everything here was produced by the pipeline in",
              "https://github.com/jinichlab/SMSTPred (`workflow/run_refresh.sh`) from the UniProt release named above; the code, the figure builders and the settings (`refresh_config.json`) are in that repository.",
              "Large files are plain `.h5` (HDF5, one vector per UniProt accession), `.csv`, `.npz` or `.tar` archives of folders. Checksums: `MANIFEST.sha256`.", "",
              "## Files", "", "| File | Size (MB) | What it is |", "|---|---:|---|"]
    readme += [f"| `{n}` | {s / 1e6:.1f} | {d} |" for n, s, d in rows]
    readme += ["", "## Notes", "",
               "- The label tables type each substrate from Rhea's own compound classification (`0_uniprot_pull/substrate_binning.py` in the repository).",
               "- The ProtT5 vectors are mean-pooled per-protein embeddings of the encoder (`Rostlab/prot_t5_xl_uniref50`); the deployed ensemble uses them (`deployment_ensemble_*.tar`).",
               "- Structures are AlphaFold DB models, or ColabFold predictions where AlphaFold DB had none or one for an older sequence; each file names its source in its header and in `structures_status.csv`.",
               "- Licence: CC BY 4.0. The label tables, sequences and structures derive from UniProt and AlphaFold DB (CC BY 4.0); ExPASy ENZYME and Rhea (CC BY 4.0) are redistributed in `reference_tables_*.tar`.",
               "- This record supersedes the earlier, unpublished draft of the same data (built from earlier UniProt releases): none of its files are kept.", ""]
    (out / "README.md").write_text("\n".join(readme))
    with open(out / "MANIFEST.sha256", "w") as f:
        for p in sorted(out.iterdir()):
            if p.is_file() and p.name != "MANIFEST.sha256":
                f.write(f"{sha256(p)}  {p.name}\n")
    print(f"\n{len(rows)} files, {sum(s for _, s, _ in rows) / 1e9:.2f} GB -> {out}")


if __name__ == "__main__":
    main()
