#!/usr/bin/env python3
"""Stage 0 -- the two public reference tables the dark-cluster annotation (Stage 6) names enzymes and reactions with.

  python fetch_reference_tables.py [--release 2026_03] [--from-dir DIR] [--force]

Writes into <release>/reference_tables/ (release = the newest in $UNIPROT_PULL_DIR unless --release is given):
  enzyme.dat    ExPASy ENZYME: every EC number with its recommended name (https://ftp.expasy.org/databases/enzyme/enzyme.dat, ~10 MB)
  rhea2ec.tsv   Rhea reaction -> EC links (https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv, ~1 MB)
  rhea-directions.tsv   Rhea master reaction -> its left-to-right, right-to-left and bidirectional entries (.../rhea/tsv/rhea-directions.tsv): the only reliable way from a directional entry to its master
  rhea_equations.tsv    every Rhea reaction's equation text (RHEA_ID, equation) from Rhea's SPARQL endpoint (https://sparql.rhea-db.org/sparql, ~75,000 rows, a few MB): the names that
                        the Fig 8 titles and the dark-cluster tables are written from
  reference_tables.json   source URL, size, sha256 and the release string found in the file (ENZYME states its own release date)
A file that is already there is kept (--force downloads it again). --from-dir copies the files from a local folder instead of the network (for an offline run, or to
seed a release with copies you already have). Network: three requests, a few seconds.
"""
import argparse
import hashlib
import json
import os
import re
import shutil
import sys
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "1_clustering"))
from mmseqs_clusters import newest_release  # noqa: E402

SPARQL = "https://sparql.rhea-db.org/sparql"
EQ_QUERY = """PREFIX rh:<http://rdf.rhea-db.org/>
SELECT ?acc ?eq WHERE { ?r rh:accession ?acc ; rh:equation ?eq . } ORDER BY ?acc LIMIT 20000 OFFSET %d"""
FILES = {"enzyme.dat": "https://ftp.expasy.org/databases/enzyme/enzyme.dat",
         "rhea2ec.tsv": "https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv",
         "rhea-directions.tsv": "https://ftp.expasy.org/databases/rhea/tsv/rhea-directions.tsv"}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def check(name, p):
    """A truncated or HTML-error download must not be taken for the table."""
    head = open(p, errors="replace").read(4000)
    if name == "enzyme.dat":
        ok = "ENZYME nomenclature database" in head
    elif name == "rhea2ec.tsv":
        ok = head.startswith("RHEA_ID\tDIRECTION\tMASTER_ID\tID")
    else:
        ok = head.startswith("RHEA_ID_MASTER\tRHEA_ID_LR\tRHEA_ID_RL\tRHEA_ID_BI")
    if not ok:
        raise SystemExit(f"{p} does not look like {name}: not kept")


def fetch_equations(dst):
    """rhea_equations.tsv (RHEA_ID<TAB>equation) by paging the SPARQL endpoint (the endpoint caps one answer)."""
    import json as _json
    rows, off = [], 0
    while True:
        req = urllib.request.Request(SPARQL + "?" + urllib.parse.urlencode({"query": EQ_QUERY % off, "format": "application/sparql-results+json"}), headers={"User-Agent": "SMSTPred-refresh"})
        with urllib.request.urlopen(req, timeout=300) as r:
            b = _json.load(r)["results"]["bindings"]
        rows += [(x["acc"]["value"].replace("RHEA:", ""), x["eq"]["value"].replace("\t", " ")) for x in b]
        if len(b) < 20000:
            break
        off += 20000
    if len(rows) < 30000:
        raise SystemExit(f"only {len(rows)} Rhea equations came back: not kept")
    tmp = dst.with_suffix(".part")
    tmp.write_text("RHEA_ID\tequation\n" + "\n".join(f"{a}\t{e}" for a, e in rows) + "\n")
    tmp.replace(dst)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", default=None)
    ap.add_argument("--from-dir", default=None, help="copy the files from this folder instead of downloading")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    pull = Path(os.environ.get("UNIPROT_PULL_DIR", HERE.parent / "data" / "uniprot_pull"))
    rel = newest_release(pull, a.release)
    out = Path(rel) / "reference_tables"
    out.mkdir(exist_ok=True)
    info = {}
    for name, url in FILES.items():
        dst = out / name
        if dst.exists() and not a.force:
            print(f"{name}: already in {out}")
        else:
            tmp = out / (name + ".part")
            if a.from_dir:
                shutil.copyfile(Path(a.from_dir) / name, tmp)
            else:
                req = urllib.request.Request(url, headers={"User-Agent": "SMSTPred-refresh"})
                with urllib.request.urlopen(req, timeout=300) as r, open(tmp, "wb") as f:
                    shutil.copyfileobj(r, f)
            check(name, tmp)
            tmp.replace(dst)
            print(f"{name}: {'copied from ' + a.from_dir if a.from_dir else 'downloaded'} ({dst.stat().st_size:,} bytes)")
        rec = {"source": url, "bytes": dst.stat().st_size, "sha256": sha256(dst)}
        if name == "enzyme.dat":
            m = re.search(r"Release of ([\w-]+)", open(dst, errors="replace").read(4000))
            rec["enzyme_release"] = m.group(1) if m else None
        info[name] = rec
    eqf = out / "rhea_equations.tsv"
    if eqf.exists() and not a.force:
        print(f"rhea_equations.tsv: already in {out}")
    elif a.from_dir and (Path(a.from_dir) / "rhea_equations.tsv").exists():
        shutil.copyfile(Path(a.from_dir) / "rhea_equations.tsv", eqf)
    else:
        fetch_equations(eqf)
        print(f"rhea_equations.tsv: downloaded ({eqf.stat().st_size:,} bytes)")
    info["rhea_equations.tsv"] = {"source": SPARQL, "bytes": eqf.stat().st_size, "sha256": sha256(eqf), "rows": sum(1 for _ in open(eqf)) - 1}
    (out / "reference_tables.json").write_text(json.dumps(info, indent=1))
    print("wrote", out / "reference_tables.json")


if __name__ == "__main__":
    main()
