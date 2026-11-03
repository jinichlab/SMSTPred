"""Find the newest trained SMSTPred classifier and make sure its files are intact.

  python get_model.py [--model-dir DIR]        prints the model folder that classify.py will use

Order of the search:
  1. --model-dir / $SMSTPRED_MODEL_DIR     a folder you point at yourself (an ensemble made by 2_training/deploy_ensemble.py)
  2. model/LATEST.json                      written by the pipeline every time a new model is trained (workflow/publish_classifier.py): it names the newest release's folder under
                                            model/, with the SHA-256 of every file. The files are checked; a missing or damaged file is fetched from the `download` URL in LATEST.json
                                            (the data deposit on Zenodo) when there is one.
So after `git pull` you always get the model of the latest release; nothing else has to be done by hand.
"""
import argparse
import hashlib
import json
import os
import tarfile
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def intact(folder, files):
    bad = [n for n, s in files.items() if not (folder / n).exists() or sha256(folder / n) != s]
    return bad


def resolve(model_dir=None, quiet=False):
    """Returns (folder, info dict). info has 'release' etc. when the model comes from LATEST.json."""
    d = model_dir or os.environ.get("SMSTPRED_MODEL_DIR")
    if d:
        return Path(d), {"release": "(user-supplied folder)"}
    latest = HERE / "model" / "LATEST.json"
    if not latest.exists():
        raise SystemExit("no model: model/LATEST.json is missing. Run `git pull`, or pass --model-dir / set SMSTPRED_MODEL_DIR to a trained ensemble folder.")
    info = json.loads(latest.read_text())
    folder = HERE / "model" / info["path"]
    bad = intact(folder, info["files"])
    if bad:
        url = info.get("download")
        if not url:
            raise SystemExit(f"model files missing or damaged in {folder}: {bad}. `git pull` / `git checkout -- ready_classifier/model` restores them.")
        if not quiet:
            print(f"fetching the model of release {info['release']} from {url} ...", flush=True)
        folder.mkdir(parents=True, exist_ok=True)
        tmp = folder / "_download.tar"
        urllib.request.urlretrieve(url, tmp)
        with tarfile.open(tmp) as t:
            for m in t.getmembers():          # only plain files of the model folder, flattened
                if m.isfile() and Path(m.name).name in info["files"]:
                    m.name = Path(m.name).name
                    t.extract(m, folder)
        tmp.unlink()
        bad = intact(folder, info["files"])
        if bad:
            raise SystemExit(f"after the download these files still do not match their checksums: {bad}")
    return folder, info


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model-dir")
    a = ap.parse_args()
    f, i = resolve(a.model_dir)
    print(f, "| release", i.get("release"))
