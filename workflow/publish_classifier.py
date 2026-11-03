#!/usr/bin/env python3
"""Make the newest trained model the one the ready classifier uses.

  python workflow/publish_classifier.py [--release <dir>] [--download-url URL]

Copies the deployed ensemble of a release (<release>/deployment: the six member weights, model_config.json, thresholds.json, the applicability-domain reference and the report) into
ready_classifier/model/<release>/ and writes ready_classifier/model/LATEST.json (release, UniProt release, labels, the honest out-of-fold metrics, and the SHA-256 of every file).
Earlier release folders under ready_classifier/model/ are removed, so the folder always holds exactly the newest model. `git pull` then gives everybody the newest model; classify.py
checks the files against LATEST.json every time it runs. --download-url (the deployment tar on the data deposit, once it is public) lets classify.py re-fetch a damaged copy itself.
It only changes the working tree: committing and pushing is a deliberate step (the chain never pushes).
"""
import argparse
import datetime
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "1_clustering"))
KEEP = ["model_config.json", "thresholds.json", "ad_reference.h5", "report.txt"]


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release")
    ap.add_argument("--download-url", default=None)
    a = ap.parse_args()
    from mmseqs_clusters import newest_release
    rel = Path(a.release or os.environ.get("RELEASE_DIR") or newest_release(Path(os.environ.get("UNIPROT_PULL_DIR", REPO / "data" / "uniprot_pull")), os.environ.get("RELEASE")))
    src = rel / "deployment"
    cfg = json.load(open(src / "model_config.json"))
    members = [f"ensemble_member_{i}.pt" for i in range(cfg["n_ensemble"])]
    for f in members + KEEP:
        if not (src / f).exists():
            sys.exit(f"{src / f} not found: the deploy stage did not finish")
    tag = rel.name.replace("release_", "")
    model = REPO / "ready_classifier" / "model"
    dst = model / tag
    for old in model.iterdir():
        if old.is_dir() and old.name != tag:
            shutil.rmtree(old)
    dst.mkdir(parents=True, exist_ok=True)
    for f in members + KEEP:
        shutil.copyfile(src / f, dst / f)
    man = json.load(open(rel / "pull_manifest.json"))
    thr = json.load(open(src / "thresholds.json"))
    info = {"release": tag, "uniprot_release": man.get("uniprot_release"), "uniprot_release_date": man.get("uniprot_release_date_text", man.get("uniprot_release_date")), "path": tag,
            "created": datetime.date.today().isoformat(), "embedding": cfg.get("embedding"), "embedding_model": "Rostlab/prot_t5_xl_half_uniref50-enc (mean over residues)",
            "labels": cfg["labels"], "n_members": cfg["n_ensemble"], "trained_on": cfg.get("trained_on"), "ad_cutoff": cfg.get("ad_cutoff"),
            "out_of_fold_metrics": thr.get("oof_metrics", {}), "download": a.download_url, "files": {f: sha256(dst / f) for f in members + KEEP}}
    (model / "LATEST.json").write_text(json.dumps(info, indent=2) + "\n")
    print(f"ready classifier model: release {tag} ({sum((dst / f).stat().st_size for f in info['files']) / 1e6:.1f} MB) -> {dst}")


if __name__ == "__main__":
    main()
