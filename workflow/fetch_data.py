#!/usr/bin/env python3
"""Put the archived Zenodo data where the scripts expect it, extract the archives, and verify checksums.

  python workflow/fetch_data.py --src zenodo:10.5281/zenodo.23655998   # download the public record
  python workflow/fetch_data.py --src /path/to/files --link             # use a local copy (e.g. the unpublished draft), symlinked
  python workflow/fetch_data.py --check                                 # only verify the checksums in data/

Files land in data/ (or $SMSTPRED_DATA). Archives are extracted once (a stamp file in data/.extracted/ records it):
  embeddings_tree_of_life_*.tar, prediction_files_mlp.tar, sam_structures.tar.gz   -> data/
  motif_analysis_intermediates.tar                                                  -> repository root (7_motif_analysis/)
  train_val_loss_curves.tar                                                         -> 2_training/
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
# archive -> directory it is extracted into (the archives carry their own top-level folder)
EXTRACT = {
    "embeddings_tree_of_life_ProstT5.tar": DATA,
    "embeddings_tree_of_life_ProtT5.tar": DATA,
    "prediction_files_mlp.tar": DATA,
    "sam_structures.tar.gz": DATA,
    "motif_analysis_intermediates.tar": REPO,
    "train_val_loss_curves.tar": REPO / "2_training",
}


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
        if f.is_file():
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


def extract():
    stamp = DATA / ".extracted"
    stamp.mkdir(exist_ok=True)
    for name, dest in EXTRACT.items():
        arch, mark = DATA / name, stamp / f"{name}.ok"
        if not arch.exists():
            print(f"  (skipping {name}: not present)")
            continue
        if mark.exists():
            continue
        print(f"  extracting {name} -> {dest.relative_to(REPO) if dest != DATA else 'data'}/")
        with tarfile.open(arch) as t:
            t.extractall(dest, filter="data") if sys.version_info >= (3, 12) else t.extractall(dest)
        mark.touch()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="", help="zenodo:<doi> or a folder holding the record's files")
    ap.add_argument("--link", action="store_true", help="symlink files from a local folder instead of copying")
    ap.add_argument("--check", action="store_true", help="only verify checksums")
    ap.add_argument("--no-extract", action="store_true")
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
        extract()
    return 0 if check() else 1


if __name__ == "__main__":
    sys.exit(main())
