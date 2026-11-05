#!/usr/bin/env python3
"""
Build the structure store for the structure-based embedders (ESM-IF, SaProt): one <Entry>.pdb per benchmark SAM protein, each checked against the CURRENT
UniProt sequence, plus a status file that says where every structure came from and which proteins still need a prediction.

  python fetch_structures.py [--pull-dir DIR] [--release 2026_03] [--reuse ../data/sam_structures] [--out <release>/structures]
  python fetch_structures.py --ingest-predictions <ColabFold output folder>      # after predict_structures.sh

For every protein with `in_benchmark` in labels_SAM.csv (<= 1,000 aa; the structure arms are not run on longer proteins), in this order:
  1. <out>/<Entry>.pdb already there and its sequence equals the current UniProt sequence -> kept
  2. a --reuse folder (the archived sam_structures) has one whose sequence matches -> copied
  3. AlphaFold DB, through its API (https://alphafold.ebi.ac.uk/api/prediction/<Entry>): the model's own file URL is used (no model version is hard-coded, so
     a new AlphaFold DB release is followed), and the model is taken only if its sequence equals the current UniProt sequence
  4. otherwise the protein goes to to_predict.fasta: no AlphaFold DB model (`no_afdb`), or one built on another sequence (`afdb_sequence_mismatch`; UniProt changed the
     sequence after the model was made, and embedding that structure would silently describe a different protein)
--ingest-predictions reads the ColabFold output (the rank_001 model of each protein), accepts a prediction only if its sequence matches, and stores it as <Entry>.pdb.

Outputs in <out>: <Entry>.pdb, structures_status.csv (Entry, status, source, model_version, plddt, model_created, n_residues, note), to_predict.fasta, structures_manifest.json.
Re-running resumes: kept structures are not fetched again. Network: one small API request per protein that is not kept or reused (4 in parallel).
"""
import argparse
import hashlib
import json
import os
import re
import shutil
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
PULL = Path(os.environ.get("UNIPROT_PULL_DIR", REPO / "data" / "uniprot_pull"))
API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M",
         "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O", "MSE": "M"}
NONSTANDARD = set("XBZUOJ")


def newest_release(pull_dir, release=None):
    if release:
        return Path(pull_dir) / f"release_{release}"
    rel = sorted(p for p in Path(pull_dir).glob("release_*") if p.is_dir() and not p.name.startswith("_"))
    if not rel:
        sys.exit(f"no release_* folder in {pull_dir}")
    return rel[-1]


def pdb_sequence(path):
    """One-letter sequence of the first chain, from its CA atoms (SEQRES is missing in ColabFold files)."""
    seq, chain, last = [], None, None
    with open(path) as f:
        for line in f:
            if line.startswith("ENDMDL"):
                break
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                c = line[21]
                chain = chain or c
                if c != chain:
                    continue
                key = (line[22:27])
                if key != last:
                    seq.append(THREE.get(line[17:20], "X"))
                    last = key
    return "".join(seq)


def same_sequence(a, b):
    """a: the sequence of a structure, b: the UniProt sequence. Equal, except that
       * a non-standard residue (X, U ...) in either sequence matches anything, and
       * a structure may LEAVE OUT the unknown residues (X) of the UniProt sequence: ColabFold does, so its model of a sequence with one X is one residue shorter
         (Q2NIH4 303 aa, model 302; O51215 286 aa, model 285). Without this such a structure looks like a different protein with every later residue shifted."""
    if len(a) == len(b) and (a == b or all(x == y or x in NONSTANDARD or y in NONSTANDARD for x, y in zip(a, b))):
        return True
    b_without_unknown = b.replace("X", "")
    return b_without_unknown != b and a == b_without_unknown


def header_source(path):
    with open(path) as f:
        head = "".join(f.readline() for _ in range(3))
    return "AFDB" if "ALPHAFOLD MONOMER" in head.upper() else "other (ColabFold?)"


def api_record(acc, session, tries=6):
    """-> (status_code, record or None). 200 with the canonical model, 404 when AlphaFold DB has none, -1 when it could not be reached (the last HTTP status goes in the note)."""
    last = "no response"
    for k in range(tries):
        try:
            r = session.get(API.format(acc=acc), timeout=60)
            if r.status_code == 200:
                recs = r.json()
                rec = next((x for x in recs if x.get("entryId") == f"AF-{acc}-F1"), recs[0] if recs else None)
                return 200, rec
            if r.status_code == 404:
                return 404, None
            last = f"HTTP {r.status_code}"
            wait = float(r.headers.get("Retry-After", 0) or 0)  # rate limit: do what the server asks
        except (requests.RequestException, ValueError) as e:
            last, wait = type(e).__name__, 0
        time.sleep(max(wait, 3 * (k + 1)))
    return -1, last


def fetch_one(args):
    acc, seq, out, session = args
    code, rec = api_record(acc, session)
    if code == 404 or (code == 200 and rec is None):
        return acc, {"status": "no_afdb", "note": "no AlphaFold DB entry"}
    if code != 200:
        return acc, {"status": "error", "note": f"AlphaFold DB API unreachable ({rec}); run again"}
    info = {"model_version": rec.get("latestVersion"), "plddt": rec.get("globalMetricValue"), "model_created": rec.get("modelCreatedDate"), "source": "AFDB"}
    if not same_sequence(rec.get("uniprotSequence", ""), seq):
        return acc, {**info, "status": "afdb_sequence_mismatch", "note": f"AlphaFold DB model is for a {len(rec.get('uniprotSequence', ''))} aa sequence, UniProt now has {len(seq)} aa"}
    last = "no response"
    for k in range(5):  # the file server throttles bursts (HTTP 403 / 429 for a while): back off and try again
        try:
            r = session.get(rec["pdbUrl"], timeout=120)
            if r.status_code == 200 and r.content.startswith((b"HEADER", b"TITLE", b"REMARK", b"MODEL", b"ATOM")):
                tmp = out / f".{acc}.pdb.tmp"
                tmp.write_bytes(r.content)
                os.replace(tmp, out / f"{acc}.pdb")
                time.sleep(0.25)  # pace the downloads
                return acc, {**info, "status": "afdb", "note": ""}
            last = f"HTTP {r.status_code}"
            wait = float(r.headers.get("Retry-After", 0) or 0)
        except (requests.RequestException, KeyError) as e:
            last, wait = type(e).__name__, 0
        time.sleep(max(wait, 10 * (k + 1)))
    return acc, {**info, "status": "error", "note": f"download of the model file failed ({last}); run again"}


def ingest(pred_dir, rows, seqs, out):
    """Take the rank_001 model of every ColabFold prediction whose sequence matches the current UniProt sequence."""
    n = 0
    for pdb in sorted(Path(pred_dir).glob("*rank_001*.pdb")):
        acc = pdb.name.split("_")[0]
        if acc not in seqs:
            continue
        if not same_sequence(pdb_sequence(pdb), seqs[acc]):
            rows[acc] = {"status": "prediction_sequence_mismatch", "note": f"{pdb.name} does not match the UniProt sequence"}
            continue
        plddt = None
        for js in Path(pred_dir).glob(f"{acc}_scores_rank_001*.json"):
            p = json.loads(js.read_text()).get("plddt")
            plddt = round(sum(p) / len(p), 2) if p else None
        shutil.copyfile(pdb, out / f"{acc}.pdb")
        rows[acc] = {"status": "colabfold", "source": "ColabFold (rank_001 by pLDDT)", "plddt": plddt, "model_version": pdb.name.split("rank_001_")[-1].replace(".pdb", ""), "note": ""}
        n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pull-dir", default=str(PULL))
    ap.add_argument("--release")
    ap.add_argument("--flag", default="in_benchmark")
    ap.add_argument("--out", help="structure store (default <release>/structures)")
    ap.add_argument("--reuse", action="append", default=[], help="folder of <Entry>.pdb to reuse where the sequence still matches (the archived sam_structures); repeatable")
    ap.add_argument("--ingest-predictions", help="ColabFold output folder to ingest, instead of fetching")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--limit", type=int, default=0, help="only the first N proteins (smoke test)")
    a = ap.parse_args()
    cache = Path(a.pull_dir) / "_structure_cache"      # every verified structure of every release (AlphaFold DB models and ColabFold predictions): reused where the sequence still matches
    if cache.is_dir() and str(cache) not in a.reuse:
        a.reuse.append(str(cache))

    rel = newest_release(a.pull_dir, a.release)
    out = Path(a.out) if a.out else rel / "structures"
    out.mkdir(parents=True, exist_ok=True)
    lab = pd.read_csv(rel / "labels_SAM.csv", usecols=["Entry", a.flag], low_memory=False)
    entries = list(lab.loc[lab[a.flag].astype(bool), "Entry"])
    if a.limit:
        entries = entries[:a.limit]
    sq = pd.read_csv(rel / "sam_mtases.tsv.gz", sep="\t", usecols=["Entry", "Sequence"], low_memory=False).drop_duplicates("Entry").set_index("Entry")["Sequence"]
    seqs = {e: str(sq[e]).strip().upper() for e in entries}
    print(f"release {rel.name}: {len(entries)} proteins with {a.flag}; structure store {out}", flush=True)

    status_file = out / "structures_status.csv"
    rows = {}
    if status_file.exists():
        rows = {r["Entry"]: {k: v for k, v in r.items() if k != "Entry" and pd.notna(v)} for r in pd.read_csv(status_file).to_dict("records")}

    if a.ingest_predictions:
        n = ingest(a.ingest_predictions, rows, seqs, out)
        print(f"ingested {n} ColabFold predictions", flush=True)
    else:
        todo, session = [], requests.Session()
        for e in entries:
            dest = out / f"{e}.pdb"
            if dest.exists() and same_sequence(pdb_sequence(dest), seqs[e]):
                rows.setdefault(e, {}).update({"status": rows.get(e, {}).get("status") if rows.get(e, {}).get("status") in ("afdb", "colabfold", "reused_archived") else "kept",
                                               "note": rows.get(e, {}).get("note", "")})
                rows[e].setdefault("source", header_source(dest))
                continue
            if dest.exists():
                dest.unlink()  # a stale file for an older sequence
            for folder in a.reuse:
                old = Path(folder) / f"{e}.pdb"
                if old.exists() and same_sequence(pdb_sequence(old), seqs[e]):
                    shutil.copyfile(old, dest)
                    rows[e] = {"status": "reused_archived", "source": header_source(dest), "note": ""}
                    break
            else:
                stale = [f for f in a.reuse if (Path(f) / f"{e}.pdb").exists()]
                todo.append((e, seqs[e], out, session))
                if stale:
                    rows[e] = {"note": "the archived structure is for an older sequence"}
        have = sum(1 for e in entries if (out / f"{e}.pdb").exists())
        print(f"already in the store (verified) {have} | reused from the archive {sum(1 for r in rows.values() if r.get('status') == 'reused_archived')} | "
              f"to ask AlphaFold DB for: {len(todo)}", flush=True)
        with ThreadPoolExecutor(a.workers) as ex:
            for acc, res in ex.map(fetch_one, todo):
                notes = [x for x in (rows.get(acc, {}).get("note", ""), res.get("note", "")) if x]
                rows[acc] = {**res, "note": "; ".join(dict.fromkeys(notes))}  # keep a note from this run once (a re-run must not repeat it)

    # the final table, and what still needs a prediction
    table = []
    for e in entries:
        r = rows.get(e, {"status": "missing", "note": "not processed"})
        has = (out / f"{e}.pdb").exists()
        n_res = len(pdb_sequence(out / f"{e}.pdb")) if has else None
        table.append({"Entry": e, "status": r.get("status", "missing") if has or r.get("status") in ("no_afdb", "afdb_sequence_mismatch", "error", "prediction_sequence_mismatch") else "missing",
                      "source": r.get("source", ""), "model_version": r.get("model_version", ""), "plddt": r.get("plddt", ""), "model_created": r.get("model_created", ""),
                      "n_residues": n_res, "note": r.get("note", "")})
    tab = pd.DataFrame(table)
    tab.to_csv(status_file, index=False)
    need = tab[~tab["Entry"].map(lambda e: (out / f"{e}.pdb").exists())]
    unreachable = need[need["status"].isin(["error", "missing"])]  # not known to need a prediction: run again
    need = need[~need["Entry"].isin(unreachable["Entry"])]
    with open(out / "to_predict.fasta", "w") as f:
        for e in need["Entry"]:
            f.write(f">{e}\n{seqs[e]}\n")
    counts = tab["status"].value_counts().to_dict()
    (out / "structures_manifest.json").write_text(json.dumps({
        "release": rel.name, "flag": a.flag, "proteins": len(entries), "status_counts": counts, "to_predict": len(need), "unreachable": len(unreachable),
        "afdb_models": sorted({str(v) for v in tab.loc[tab["status"] == "afdb", "model_version"]}), "reuse": a.reuse,
        "date": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "api": "alphafold.ebi.ac.uk/api/prediction (the file URL of each record is used; no model version is hard-coded)"}, indent=2))
    print(f"status: {counts}")
    print(f"{len(entries) - len(need) - len(unreachable)} of {len(entries)} proteins have a verified structure; {len(need)} need a prediction -> {out / 'to_predict.fasta'}"
          + (f"; {len(unreachable)} could not be checked (AlphaFold DB unreachable): run again" if len(unreachable) else ""))
    if len(need):
        print("run predict_structures.sh next, then `fetch_structures.py --ingest-predictions <its output folder>`")
    cache.mkdir(exist_ok=True)         # keep what this release verified (never over a file that is there: a changed sequence is caught when a file is reused)
    new = 0
    for f in out.glob("*.pdb"):
        if not (cache / f.name).exists():
            shutil.copyfile(f, cache / f.name)
            new += 1
    print(f"structure cache {cache}: {new} new structures added, {sum(1 for _ in cache.glob('*.pdb'))} in all")
    if unreachable:       # an incomplete set must not pass as done: the next stages (and the benchmark) would silently lose these proteins; the chain runs this script again
        print(f"INCOMPLETE: {len(unreachable)} proteins could not be checked or downloaded; exiting with an error so the stage is repeated", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
