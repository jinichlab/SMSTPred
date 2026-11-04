#!/usr/bin/env python3
"""Put the archived Zenodo data where the scripts expect it, rebuild the release folder from it, and verify checksums.

  python workflow/fetch_data.py --src zenodo:10.5281/zenodo.23655998   # download the public record
  python workflow/fetch_data.py --src /path/to/files --link             # use a local copy (e.g. the unpublished draft), symlinked
  python workflow/fetch_data.py --check                                 # only verify the checksums in data/

The record's files land in data/ (or $SMSTPRED_DATA). They are then assembled into the release folder the pipeline stages read,
data/uniprot_pull/release_<tag>/ (or $UNIPROT_PULL_DIR/release_<tag>): single files are linked to their place (embeddings/, labels_SAM.csv,
predictions_all_proteomes.csv.gz, ...), the .tar archives (deployment, benchmark, umap, tree, dark_annotation, motif_analysis,
reference_tables) are extracted into it, and the structure store is extracted to structures/. This is the exact inverse of
workflow/build_zenodo_files.py (it reads the same table), so `RELEASE_DIR=<that folder>` runs any stage on the archived data.
The large regenerable motif inputs (corpus, class FASTAs) are not archived; stage 7 steps that need them are tier C.
"""
import argparse
import hashlib
import json
import os
import shutil
import sys
import tarfile
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("SMSTPRED_DATA", REPO / "data"))
sys.path.insert(0, str(REPO / "workflow"))
PULL = Path(os.environ.get("UNIPROT_PULL_DIR", DATA / "uniprot_pull"))


def sha256(path, block=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(block):
            h.update(chunk)
    return h.hexdigest()


def manifest(directory):
    m = {}
    p = Path(directory) / "MANIFEST.sha256"
    if p.exists():
        for line in p.read_text().splitlines():
            if line.strip():
                digest, name = line.split(None, 1)
                m[name.strip().lstrip("*")] = digest
    return m


def check(directory=None):
    d = Path(directory or DATA)
    m = manifest(d)
    if not m:
        print(f"no MANIFEST.sha256 in {d}")
        return False
    bad = missing = ok = 0
    for name, digest in sorted(m.items()):
        f = d / name
        if name == "README.md":      # the record's README is not placed (data/README.md is the repository's)
            continue
        if not f.exists():
            missing += 1
            print(f"  MISSING  {name}")
        elif sha256(f) != digest:
            bad += 1
            print(f"  BAD      {name}")
        else:
            ok += 1
            print(f"  OK       {name}")
    print(f"{ok} ok, {bad} bad, {missing} missing")
    return bad == 0 and missing == 0


def place_local(src, link):
    DATA.mkdir(parents=True, exist_ok=True)
    for f in sorted(Path(src).iterdir()):
        if f.is_file() and f.name != "README.md":      # the record's README would replace data/README.md
            dest = DATA / f.name
            if dest.exists() or dest.is_symlink():
                if dest.stat().st_size == f.stat().st_size:
                    continue
                dest.unlink()
            if link:
                os.symlink(f.resolve(), dest)
            else:
                shutil.copy2(f, dest)
            print(f"  {'linked' if link else 'copied'} {f.name}")


def download(doi):
    rec = doi.split("zenodo.")[-1]
    meta = json.load(urllib.request.urlopen(f"https://zenodo.org/api/records/{rec}"))
    DATA.mkdir(parents=True, exist_ok=True)
    for f in meta["files"]:
        dest = DATA / f["key"]
        if dest.exists() and dest.stat().st_size == f["size"]:
            continue
        print(f"  downloading {f['key']} ({f['size'] / 1e6:.0f} MB)")
        urllib.request.urlretrieve(f["links"]["self"], dest)


def split_proteome_embeddings(rel):
    """embeddings/embeddings_proteomes_ProtT5.h5 (one file, keyed by accession) -> embeddings/proteomes_ProtT5/<proteome>_embeddings.h5 (the per-proteome files Stage 3 reads)."""
    import h5py
    import pandas as pd
    big, outd = rel / "embeddings" / "embeddings_proteomes_ProtT5.h5", rel / "embeddings" / "proteomes_ProtT5"
    members = rel / "embedding_inputs" / "proteome_members.csv"
    if not big.exists() or not members.exists() or outd.exists():
        return
    m = pd.read_csv(members, usecols=["proteome", "Entry"])
    outd.mkdir(parents=True)
    with h5py.File(big, "r") as f:
        for prot, g in m.groupby("proteome"):
            with h5py.File(outd / f"{prot}_embeddings.h5", "w") as o:
                for e in g.Entry:
                    if e in f:
                        o.create_dataset(e, data=f[e][()])
    print(f"  split the proteome embeddings into {m.proteome.nunique()} files")


def assemble(tag=None, link=True):
    """Rebuild <UNIPROT_PULL_DIR>/release_<tag>/ from the files in DATA (inverse of build_zenodo_files.FILES)."""
    from build_zenodo_files import FILES
    if tag is None:
        names = [p.name for p in DATA.glob("labels_SAM_*.csv")]
        if not names:
            print("  no labels_SAM_<tag>.csv in data/: nothing to assemble")
            return None
        tag = names[0][len("labels_SAM_"):-len(".csv")]
    rel = PULL / f"release_{tag}"
    stamp = rel / ".assembled"
    if stamp.exists():
        return rel
    rel.mkdir(parents=True, exist_ok=True)
    for name, kind, src, _ in FILES:
        name = name.replace("RELEASE", tag)
        f = DATA / name
        if not f.exists():
            print(f"  (skipping {name}: not present)")
            continue
        if kind in ("copy", "annot"):
            dst = rel / src
            dst.parent.mkdir(parents=True, exist_ok=True)
            if dst.exists() or dst.is_symlink():
                dst.unlink()
            os.symlink(f.resolve(), dst) if link else shutil.copyfile(f, dst)
        else:
            print(f"  extracting {name} -> {rel.name}/")
            with tarfile.open(f) as t:
                t.extractall(rel, filter="data") if sys.version_info >= (3, 12) else t.extractall(rel)
    split_proteome_embeddings(rel)
    stamp.touch()
    print(f"  release folder: {rel}")
    return rel


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="", help="zenodo:<doi> or a folder holding the record's files")
    ap.add_argument("--link", action="store_true", help="symlink files from a local folder instead of copying")
    ap.add_argument("--check", action="store_true", help="only verify checksums")
    ap.add_argument("--no-extract", action="store_true", help="do not assemble the release folder")
    a = ap.parse_args()
    if a.check:
        return 0 if check() else 1
    if not a.src:
        print("give --src zenodo:10.5281/zenodo.23655998 or --src <folder with the record's files>")
        return 2
    if a.src.startswith("zenodo:"):
        download(a.src[len("zenodo:"):])
    else:
        place_local(a.src, a.link)
    if not a.no_extract:
        assemble()
    return 0 if check() else 1


if __name__ == "__main__":
    sys.exit(main())
