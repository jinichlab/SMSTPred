#!/usr/bin/env python3
"""
Stage 0 -- pull the current UniProt data the analysis starts from, through the UniProt REST API.

  python pull_uniprot.py sam         Swiss-Prot methyltransferases of the SAM-dependent MTase superfamily (InterPro IPR029063)
  python pull_uniprot.py up          Swiss-Prot enzymes with a catalytic-activity annotation (the broad "UP" comparison set)
  python pull_uniprot.py proteomes   the 1,259 reference proteomes of the tree-of-life analysis (data/selected_1259_organisms.csv)
  python pull_uniprot.py all         the three above, in that order
  python pull_uniprot.py resplit     re-split an existing pull at another --length-view (no network)

What is pulled (queries as sent to the API; they reproduce how the archived sets were defined, checked against the archived tables
-- see README.md):
  sam        (reviewed:true) AND (xref:interpro-IPR029063)
  up         (reviewed:true) AND (cc_catalytic_activity:*)            the archived UP table is 100 % reviewed, 100 % with catalytic-activity text
  proteomes  (proteome:UPxxxxxxxxx)  for every Proteome Id in selected_1259_organisms.csv

NOTHING IS FILTERED BY LENGTH at pull time. Proteins longer than --length-view (default 2,000 aa: the length the deployed model is
trained and applied on in the refresh; the ARCHIVED training sets and embeddings stop at 1,000 aa) are kept in the raw tables and
listed separately:
  sam, up      the TSV has every length; the FASTA files hold the <= 2,000 aa proteins (what gets embedded)
  proteomes    reference_proteomes/<organism>.csv      proteins <= 2,000 aa  (same file name and columns as the archived tables,
                                                       so 4_phylo_trends / 7_motif_analysis can read them unchanged)
               reference_proteomes_long/<organism>.csv proteins  > 2,000 aa  (kept, not used downstream)
  proteome_summary.csv has the counts of both.
The cutoff can be changed afterwards without downloading again: `python pull_uniprot.py resplit --length-view L`.
Methods that cannot take long proteins (SaProt: 1,024 tokens; ESM-IF and Foldseek: AlphaFold DB has no full-length model above
about 2,700 aa) are evaluated on the <= 1,000 aa subset: filter on the Length column.

Nothing already in the repository or in the archived data is overwritten: output goes to <out>/release_<UniProt release>/ and
the release is pinned (the pull stops if the release changes while it runs; --expect-release makes that explicit).

  --out DIR            default $UNIPROT_PULL_DIR, else <repo>/data/uniprot_pull
  --limit N            smoke test: N rows (sam, up) / N proteomes (proteomes), written under <out>/_smoke/
  --workers N          parallel proteome downloads (default 4; the client also spaces requests >= 0.2 s apart)
  --redo               download proteomes again that the status file already lists as done
  --expect-release R   stop unless UniProt is serving release R (e.g. 2026_03)
  --length-view L      the length cutoff for the "<= L aa" views (default 2000; 0 = no split)
  --release R          resplit only: which release_R folder to re-split (default: the newest)
  --method stream|paged  sam / up: one streamed request with a row-count check and retries (default), or cursor pages of 500 rows
"""
import argparse
import datetime as dt
import io
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from uniprot_api import Client, sha256

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
SAM_QUERY = "(reviewed:true) AND (xref:interpro-IPR029063)"
UP_QUERY = "(reviewed:true) AND (cc_catalytic_activity:*)"
# The columns of the archived SAM / UP tables, plus organism, EC, Pfam and the entry / sequence dates (to see which sequences changed).
FIELDS_SET = ("accession,reviewed,id,protein_name,organism_name,organism_id,length,sequence,annotation_score,xref_interpro,xref_pfam,"
              "cc_function,cc_catalytic_activity,rhea,ec,ft_dna_bind,xref_chembl,xref_string,date_modified,date_sequence_modified,"
              "version,sequence_version,xref_proteomes")
# The columns of the archived per-proteome tables (uniprot_repull_1259_reviewed.py).
FIELDS_PROT = "accession,reviewed,protein_name,organism_name,length,sequence,xref_interpro,xref_pfam"
KEEP = ["Entry", "Reviewed", "Protein names", "Organism", "Length", "Sequence", "InterPro", "Pfam"]
SAM_DOMAIN = "IPR029063"


def canon(org):
    """File stem of an organism, as in the archived tables and prediction files."""
    return str(org).replace(" ", "_").replace("/", "_")


def now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def say(msg):
    print(msg, flush=True)


# ---------------------------------------------------------------------------------------------- release / output folder
def setup(client, args):
    n, release, rdate = client.count(SAM_QUERY)
    if args.expect_release and release != args.expect_release:
        sys.exit(f"UniProt is serving release {release}, not the expected {args.expect_release}")
    base = Path(args.out or os.environ.get("UNIPROT_PULL_DIR") or REPO / "data" / "uniprot_pull")
    if args.limit:
        base = base / "_smoke"
    out = base / f"release_{release}"
    out.mkdir(parents=True, exist_ok=True)
    say(f"UniProt release {release} ({rdate}); writing to {out}")
    return out, release, rdate


def update_manifest(out, release, rdate, key, info, client, args):
    p = out / "pull_manifest.json"
    m = json.loads(p.read_text()) if p.exists() else {}
    m.update({"uniprot_release": release, "uniprot_release_date": rdate, "api": "https://rest.uniprot.org",
              "length_view": args.length_view})
    m.setdefault("datasets", {})[key] = {**info, "pulled_at": now(), "command": " ".join(sys.argv)}
    seen = sorted({r for r, _ in client.releases})
    m["releases_seen_in_this_session"] = seen
    if len(seen) > 1:
        m["WARNING"] = f"more than one UniProt release was served during the pull: {seen}"
        say(f"  WARNING: releases seen: {seen}")
    try:
        m["repo_commit"] = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True, text=True).stdout.strip() or None
    except Exception:
        m["repo_commit"] = None
    p.write_text(json.dumps(m, indent=1))
    if not args.limit:                                   # a small copy next to the code, to track in git
        (HERE / "pull_manifests").mkdir(exist_ok=True)
        (HERE / "pull_manifests" / f"pull_manifest_{release}.json").write_text(json.dumps(m, indent=1))


# ---------------------------------------------------------------------------------------------- sam / up
def write_fasta(tsv_gz, fasta, max_len):
    d = pd.read_csv(tsv_gz, sep="\t", usecols=["Entry", "Length", "Sequence"], dtype={"Entry": str, "Sequence": str})
    n_all = len(d)
    if max_len:
        d = d[d["Length"] <= max_len]
    with open(fasta, "w") as f:
        for e, s in zip(d["Entry"], d["Sequence"]):
            f.write(f">{e}\n{s}\n")
    return n_all, len(d)


def pull_set(client, key, query, out, release, rdate, args):
    say(f"\n== {key}: {query}")
    n, _, _ = client.count(query)
    say(f"  {n} entries (all lengths)")
    tsv = out / f"{key}.tsv.gz"
    if args.limit or args.method == "paged":
        res = client.paged_tsv(query, FIELDS_SET, tsv, limit=args.limit, log=say)       # a --limit smoke test stops after a few pages
    else:
        res = client.stream_tsv_gz(query, FIELDS_SET, tsv, log=say)
    fasta = out / f"{key}_le{args.length_view}.fasta" if args.length_view else out / f"{key}.fasta"
    n_all, n_fa = write_fasta(tsv, fasta, args.length_view)
    say(f"  {res['rows']} rows in {res['seconds']} s -> {tsv.name}; FASTA {fasta.name}: {n_fa} proteins <= {args.length_view} aa")
    info = {"query": query, "fields": FIELDS_SET, "rows": res["rows"], "rows_announced_by_api": res["expected"],
            "rows_le_length_view": n_fa, "rows_gt_length_view": n_all - n_fa, "seconds": res["seconds"],
            "release_of_download": res["release"],
            "files": {tsv.name: {"bytes": tsv.stat().st_size, "sha256": sha256(tsv)},
                      fasta.name: {"bytes": fasta.stat().st_size, "sha256": sha256(fasta)}}}
    update_manifest(out, release, rdate, key, info, client, args)
    return info


# ---------------------------------------------------------------------------------------------- proteomes
def fetch_one(client, pid, organism, stem, out, length_view):
    t0 = time.time()
    rec = {"proteome_id": pid, "organism": organism, "stem": stem}
    try:
        rec.update(client.proteome_meta(pid))
        text, rel = client.stream_tsv(f"(proteome:{pid})", FIELDS_PROT)
        rec["release"] = rel
        if text.count("\n") < 2:
            rec.update(status="empty", n_all=0)
            return rec
        d = pd.read_csv(io.StringIO(text), sep="\t", low_memory=False)
        missing = [c for c in KEEP if c not in d.columns]
        if missing:
            rec.update(status=f"error: columns missing {missing}")
            return rec
        d = d[KEEP]
        ipr = d["InterPro"].fillna("").str.contains(SAM_DOMAIN, regex=False)
        short = (d["Length"] <= length_view) if length_view else pd.Series(True, index=d.index)
        for name, sub in (("reference_proteomes", d[short]), ("reference_proteomes_long", d[~short])):
            folder = out / name
            folder.mkdir(exist_ok=True)
            tmp, final = folder / f".{stem}.part", folder / f"{stem}.csv"
            if len(sub) or name == "reference_proteomes":
                sub.to_csv(tmp, index=False)
                os.replace(tmp, final)
            elif final.exists():
                final.unlink()                        # a re-pull with no long proteins left
        rec.update(status="ok", n_all=len(d), n_le=int(short.sum()), n_gt=int((~short).sum()),
                   n_reviewed_all=int((d["Reviewed"] == "reviewed").sum()), n_ipr_all=int(ipr.sum()),
                   n_ipr_le=int((ipr & short).sum()), n_ipr_gt=int((ipr & ~short).sum()))
    except Exception as e:                            # recorded, and the run goes on with the next proteome
        rec.update(status=f"error: {type(e).__name__}: {str(e)[:200]}")
    rec["seconds"] = round(time.time() - t0, 1)
    return rec


def pull_proteomes(client, out, release, rdate, args):
    org_csv = Path(args.organisms or REPO / "data" / "selected_1259_organisms.csv")
    sel = pd.read_csv(org_csv)
    if args.limit:
        sel = sel.head(args.limit)
    stems = sel["Organism"].map(canon)
    if stems.duplicated().any():                      # two proteomes with the same organism name: keep the Proteome Id in the name
        stems = stems.where(~stems.duplicated(keep=False), stems + "__" + sel["Proteome Id"])
    say(f"\n== proteomes: {len(sel)} from {org_csv.name} -> {out / 'reference_proteomes'}")
    status = out / "proteome_status.jsonl"
    done = {}
    if status.exists() and not args.redo:
        for line in status.read_text().splitlines():
            r = json.loads(line)
            if r["status"] in ("ok", "empty"):
                done[r["proteome_id"]] = r
        say(f"  resuming: {len(done)} proteomes already done")
    todo = [(p, o, s) for p, o, s in zip(sel["Proteome Id"], sel["Organism"], stems) if p not in done]
    t0, n_done, n_err = time.time(), 0, 0
    with open(status, "a") as log, ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(fetch_one, client, p, o, s, out, args.length_view) for p, o, s in todo]
        for f in as_completed(futs):
            r = f.result()
            log.write(json.dumps(r) + "\n")
            log.flush()
            if r["status"] in ("ok", "empty"):
                done[r["proteome_id"]] = r
            else:
                n_err += 1
                say(f"  [{r['status'][:90]}] {r['stem']}")
            n_done += 1
            if n_done % 50 == 0 or n_done == len(todo):
                say(f"  {n_done}/{len(todo)} downloaded ({n_err} errors), {time.time() - t0:.0f} s")
    # summary of the latest record per proteome, in the order of the organism table
    latest = {}
    for line in status.read_text().splitlines():
        r = json.loads(line)
        latest[r["proteome_id"]] = r
    rows = [latest[p] for p in sel["Proteome Id"] if p in latest]
    summ = pd.DataFrame(rows)
    summ.to_csv(out / "proteome_summary.csv", index=False)
    ok = summ[summ.status == "ok"]
    say(f"\n  ok {len(ok)} | empty {int((summ.status == 'empty').sum())} | error {int(summ.status.str.startswith('error').sum())} "
        f"of {len(sel)}")
    say(f"  proteins: {int(ok.n_all.sum())} (<= {args.length_view} aa: {int(ok.n_le.sum())}, longer: {int(ok.n_gt.sum())}); "
        f"IPR029063: {int(ok.n_ipr_all.sum())} (<= {args.length_view} aa: {int(ok.n_ipr_le.sum())})")
    empties = summ[summ.status == "empty"]
    if len(empties):
        say("  proteomes with no UniProtKB entries now (see proteome_type / proteome_redundant_to in proteome_summary.csv):")
        for _, r in empties.head(15).iterrows():
            say(f"    {r.proteome_id} {r.organism}: type={r.get('proteome_type')} redundant_to={r.get('proteome_redundant_to')}")
    info = {"query": "(proteome:<Proteome Id>) for every Proteome Id of the organism table", "fields": FIELDS_PROT,
            "organism_table": org_csv.name, "proteomes_requested": len(sel), "proteomes_ok": int(len(ok)),
            "proteomes_empty": int((summ.status == "empty").sum()),
            "proteomes_error": int(summ.status.str.startswith("error").sum()),
            "proteins_all": int(ok.n_all.sum()), "proteins_le_length_view": int(ok.n_le.sum()),
            "proteins_gt_length_view": int(ok.n_gt.sum()), "ipr029063_all": int(ok.n_ipr_all.sum()),
            "ipr029063_le_length_view": int(ok.n_ipr_le.sum()),
            "releases_of_downloads": sorted({str(x) for x in summ.get("release", pd.Series(dtype=str)).dropna()}),
            "files": {"proteome_summary.csv": {"bytes": (out / "proteome_summary.csv").stat().st_size,
                                               "sha256": sha256(out / "proteome_summary.csv")}}}
    update_manifest(out, release, rdate, "proteomes", info, client, args)
    if info["proteomes_error"]:
        say(f"  {info['proteomes_error']} proteomes failed: run again (it resumes) before using the data")
    return info


# ---------------------------------------------------------------------------------------------- main
def _resplit_one(job):
    """Re-split one organism: the short and the long file together are the whole proteome, so nothing is lost or downloaded again."""
    stem, out, view = job
    short_f, long_f = out / "reference_proteomes" / f"{stem}.csv", out / "reference_proteomes_long" / f"{stem}.csv"
    d = pd.concat([pd.read_csv(f, low_memory=False) for f in (short_f, long_f) if f.exists()], ignore_index=True)
    ipr = d["InterPro"].fillna("").str.contains(SAM_DOMAIN, regex=False)
    short = (d["Length"] <= view) if view else pd.Series(True, index=d.index)
    for name, sub in (("reference_proteomes", d[short]), ("reference_proteomes_long", d[~short])):
        folder = out / name
        folder.mkdir(exist_ok=True)
        tmp, final = folder / f".{stem}.part", folder / f"{stem}.csv"
        if len(sub) or name == "reference_proteomes":
            sub.to_csv(tmp, index=False)
            os.replace(tmp, final)
        elif final.exists():
            final.unlink()
    return {"stem": stem, "n_all": len(d), "n_le": int(short.sum()), "n_gt": int((~short).sum()),
            "n_ipr_all": int(ipr.sum()), "n_ipr_le": int((ipr & short).sum()), "n_ipr_gt": int((ipr & ~short).sum())}


def resplit(out, args):
    """Change the length cutoff of an existing pull (FASTA files, per-proteome tables, summary, manifest); no network."""
    mpath = out / "pull_manifest.json"
    m = json.loads(mpath.read_text())
    old, view = m.get("length_view"), args.length_view
    say(f"re-splitting {out}: length view {old} -> {view}")
    fname = lambda key, v: f"{key}_le{v}.fasta" if v else f"{key}.fasta"
    for key in ("sam_mtases", "up_enzymes"):
        tsv = out / f"{key}.tsv.gz"
        if not tsv.exists():
            continue
        n_all, n_fa = write_fasta(tsv, out / fname(key, view), view)
        if fname(key, old) != fname(key, view) and (out / fname(key, old)).exists():
            (out / fname(key, old)).unlink()
        ds = m["datasets"][key]
        ds["rows_le_length_view"], ds["rows_gt_length_view"] = n_fa, n_all - n_fa
        ds["files"].pop(fname(key, old), None)
        ds["files"][fname(key, view)] = {"bytes": (out / fname(key, view)).stat().st_size, "sha256": sha256(out / fname(key, view))}
        say(f"  {key}: {n_fa} of {n_all} proteins <= {view} aa -> {fname(key, view)}")
    status = out / "proteome_status.jsonl"
    if status.exists():
        latest = {}
        for line in status.read_text().splitlines():
            r = json.loads(line)
            latest[r["proteome_id"]] = r
        todo = [r for r in latest.values() if r["status"] == "ok"]
        t0, done = time.time(), 0
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            for res in ex.map(_resplit_one, [(r["stem"], out, view) for r in todo], chunksize=8):
                latest_by_stem = next(r for r in todo if r["stem"] == res["stem"])
                latest_by_stem.update({k: v for k, v in res.items() if k != "stem"})
                done += 1
                if done % 200 == 0:
                    say(f"  {done}/{len(todo)} proteomes re-split, {time.time() - t0:.0f} s")
        status.write_text("".join(json.dumps(r) + "\n" for r in latest.values()))
        sel = pd.read_csv(Path(args.organisms or REPO / "data" / "selected_1259_organisms.csv"))
        summ = pd.DataFrame([latest[p] for p in sel["Proteome Id"] if p in latest])
        summ.to_csv(out / "proteome_summary.csv", index=False)
        ok = summ[summ.status == "ok"]
        say(f"  proteins <= {view} aa: {int(ok.n_le.sum())} of {int(ok.n_all.sum())}; IPR029063 <= {view} aa: {int(ok.n_ipr_le.sum())} of {int(ok.n_ipr_all.sum())}")
        pi = m["datasets"]["proteomes"]
        pi.update({"proteins_le_length_view": int(ok.n_le.sum()), "proteins_gt_length_view": int(ok.n_gt.sum()),
                   "ipr029063_le_length_view": int(ok.n_ipr_le.sum())})
        pi["files"]["proteome_summary.csv"] = {"bytes": (out / "proteome_summary.csv").stat().st_size,
                                               "sha256": sha256(out / "proteome_summary.csv")}
    m["length_view"] = view
    m.setdefault("resplit_history", []).append({"from": old, "to": view, "at": now()})
    mpath.write_text(json.dumps(m, indent=1))
    if "_smoke" not in str(out):
        (HERE / "pull_manifests").mkdir(exist_ok=True)
        (HERE / "pull_manifests" / f"pull_manifest_{m['uniprot_release']}.json").write_text(json.dumps(m, indent=1))
    say("done.")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=["sam", "up", "proteomes", "all", "resplit"])
    ap.add_argument("--release", help="resplit only: the release_<R> folder to re-split (default: the newest)")
    ap.add_argument("--out")
    ap.add_argument("--organisms")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--redo", action="store_true")
    ap.add_argument("--expect-release")
    ap.add_argument("--length-view", type=int, default=2000)
    ap.add_argument("--method", choices=["stream", "paged"], default="stream",
                    help="how sam / up are downloaded: one streamed request (default) or 500-row pages (slower, resumable by hand)")
    args = ap.parse_args()

    if args.what == "resplit":                         # no network: just the existing release folder
        base = Path(args.out or os.environ.get("UNIPROT_PULL_DIR") or REPO / "data" / "uniprot_pull")
        if args.limit:
            base = base / "_smoke"
        rels = sorted(p for p in base.glob("release_*") if p.is_dir())
        if not rels:
            sys.exit(f"no release_* folder in {base}")
        resplit(base / f"release_{args.release}" if args.release else rels[-1], args)
        return

    client = Client()
    out, release, rdate = setup(client, args)
    if args.what in ("sam", "all"):
        pull_set(client, "sam_mtases", SAM_QUERY, out, release, rdate, args)
    if args.what in ("up", "all"):
        pull_set(client, "up_enzymes", UP_QUERY, out, release, rdate, args)
    status_ok = True
    if args.what in ("proteomes", "all"):
        info = pull_proteomes(client, out, release, rdate, args)
        status_ok = info["proteomes_error"] == 0
    seen = sorted({r for r, _ in client.releases})
    say(f"\nUniProt releases seen: {seen}")
    if len(seen) > 1:
        sys.exit("the UniProt release changed during the pull: delete the release folder you do not want and pull again")
    if not status_ok:
        sys.exit(2)
    say(f"done. output: {out}")


if __name__ == "__main__":
    main()
