#!/usr/bin/env python3
"""
Due-diligence for the apo/AlphaFold-only motifs.

For the motifs where the representative had no substrate/cofactor-bound structure
(so we settled for an apo crystal or an AlphaFold model), search the FULL clusters:
every cluster MEMBER carrying the motif (>=50% id to the rep, so same fold + motif)
is checked for a bound PDB structure we could use instead.

Pipeline
  1. per apo/AF motif: FIMO-matching reps (fimo_<class>/fimo.tsv) -> expand to all
     cluster members (clust_cluster.tsv, matched on bare accession).
  2. UniProt accession -> PDB   (UniProt id-mapping REST API).
  3. each PDB -> ligand / polymer composition (RCSB data API) -> classify
     apo | cofactor-bound (SAM/SAH/analog) | substrate-bound (nucleic / other ligand).
  4. report, per motif, whether a bound structure exists among the cluster members.

Requires outbound HTTPS to rest.uniprot.org and data.rcsb.org.
Outputs -> structure_search_output/.
"""
import csv
import json
import time
import re
import urllib.request
import urllib.parse
import collections
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "structure_search_output"
OUT.mkdir(exist_ok=True)
CLUST = HERE / "clust_cluster.tsv"

# (class, FIMO motif number): current apo/AF placement being second-guessed
APO = {
    ("RNA", "1"): "m5C catalytic PC  (now: apo RsmF 3M6U)",
    ("RNA", "4"): "novel             (now: apo RsmF 3M6U)",
    ("RNA", "5"): "novel             (now: apo RsmF 3M6U)",
    ("RNA", "8"): "novel             (now: AlphaFold RlmE)",
    ("Protein", "2"): "Motif I/post-I (now: AlphaFold)",
    ("Protein", "3"): "Motif II       (now: AlphaFold)",
    ("Protein", "5"): "THW loop       (now: AlphaFold)",
    ("Protein", "8"): "NPPY catalytic (now: AlphaFold)",
}

COFACTORS = {"SAM", "SAH", "SFG", "MTA", "DZF", "SX0", "EEM", "5AD", "SAE", "SMM"}
IGNORE = {"HOH", "GOL", "EDO", "PEG", "PG4", "PG0", "1PE", "P6G", "SO4", "PO4",
          "CL", "NA", "MG", "K", "ZN", "CA", "MN", "NI", "FE", "CO", "CD", "CU",
          "ACT", "DMS", "IOD", "BR", "FMT", "MPD", "TRS", "EPE", "BME", "ACE",
          "NO3", "CAC", "IMD", "BCT", "MES", "FLC", "GG"}


def build_candidates():
    clu = collections.defaultdict(list)
    for line in open(CLUST):
        rep, mem = line.rstrip("\n").split("\t")
        clu[rep.split("|")[0]].append(mem.split("|")[0])
    per_motif, allacc = {}, set()
    for (cls, mnum), _desc in APO.items():
        reps = set()
        with open(HERE / f"fimo_{cls}" / "fimo.tsv") as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if not row or row[0].startswith("#") or row[0] == "motif_id":
                    continue
                if row[0].split("-")[0] == mnum:
                    reps.add(row[2].split("|")[0])
        accs = set(reps)
        for r in reps:
            accs.update(clu.get(r, []))
        per_motif[(cls, mnum)] = accs
        allacc |= accs
    return per_motif, sorted(allacc)


def uniprot_to_pdb(accs):
    data = urllib.parse.urlencode({"from": "UniProtKB_AC-ID", "to": "PDB",
                                   "ids": ",".join(accs)}).encode()
    job = json.load(urllib.request.urlopen(
        urllib.request.Request("https://rest.uniprot.org/idmapping/run", data=data)))["jobId"]
    print(f"[uniprot] job {job}", flush=True)
    for _ in range(200):
        r = json.load(urllib.request.urlopen(
            f"https://rest.uniprot.org/idmapping/status/{job}"))
        if r.get("jobStatus") == "FINISHED" or "results" in r:
            break
        time.sleep(3)
    url = f"https://rest.uniprot.org/idmapping/results/{job}?format=tsv&size=500"
    acc2pdb = collections.defaultdict(set)
    while url:
        resp = urllib.request.urlopen(url)
        for line in resp.read().decode().splitlines()[1:]:
            a, p = line.split("\t")[:2]
            acc2pdb[a].add(p)
        m = re.search(r'<([^>]+)>;\s*rel="next"', resp.headers.get("Link", ""))
        url = m.group(1) if m else None
    return acc2pdb


def rcsb_entry(pdb, cache={}):
    if pdb in cache:
        return cache[pdb]
    for attempt in range(3):
        try:
            d = json.load(urllib.request.urlopen(
                f"https://data.rcsb.org/rest/v1/core/entry/{pdb}", timeout=30))
            break
        except Exception:
            time.sleep(1.5)
    else:
        cache[pdb] = None
        return None
    info = d.get("rcsb_entry_info", {})
    ligs = [x for x in (info.get("nonpolymer_bound_components") or [])]
    poly = (info.get("polymer_composition") or "")
    has_na = "NA" in poly or "nucleic" in poly.lower()
    cof = [l for l in ligs if l in COFACTORS]
    subs = [l for l in ligs if l not in IGNORE and l not in COFACTORS]
    tags = []
    if has_na:
        tags.append("nucleic-substrate")
    if cof:
        tags.append("cofactor:" + "+".join(cof))
    if subs:
        tags.append("ligand:" + "+".join(subs[:4]))
    rec = {"pdb": pdb, "title": d.get("struct", {}).get("title", "")[:90],
           "resolution": (info.get("resolution_combined") or [None])[0],
           "polymer": poly, "status": "bound" if tags else "apo",
           "detail": "; ".join(tags)}
    cache[pdb] = rec
    time.sleep(0.12)
    return rec


def main():
    per_motif, allacc = build_candidates()
    print(f"[candidates] {len(allacc):,} unique cluster-member accessions", flush=True)
    with open(OUT / "candidate_accessions.csv", "w") as fh:
        fh.write("class,motif,n_accessions\n")
        for (cls, mnum), accs in per_motif.items():
            fh.write(f"{cls},{mnum},{len(accs)}\n")

    acc2pdb = uniprot_to_pdb(allacc)
    with open(OUT / "uniprot_pdb_map.tsv", "w") as fh:
        fh.write("accession\tpdb\n")
        for a, ps in acc2pdb.items():
            for p in ps:
                fh.write(f"{a}\t{p}\n")
    pdbs = sorted({p for ps in acc2pdb.values() for p in ps})
    print(f"[uniprot] {len(acc2pdb)} accessions have a PDB; {len(pdbs)} unique PDBs", flush=True)

    entries = {}
    for i, p in enumerate(pdbs, 1):
        entries[p] = rcsb_entry(p)
        if i % 25 == 0:
            print(f"[rcsb] {i}/{len(pdbs)}", flush=True)
    with open(OUT / "pdb_ligands.csv", "w") as fh:
        w = csv.writer(fh)
        w.writerow(["pdb", "status", "detail", "resolution", "polymer", "title"])
        for p in pdbs:
            e = entries[p]
            if e:
                w.writerow([e["pdb"], e["status"], e["detail"], e["resolution"],
                            e["polymer"], e["title"]])

    # per-motif report
    lines = ["Bound-structure search for apo/AF-only motifs",
             "=" * 60]
    hit_rows = []
    for (cls, mnum), desc in APO.items():
        accs = per_motif[(cls, mnum)]
        hits = []
        for a in accs:
            for p in acc2pdb.get(a, []):
                e = entries.get(p)
                if e:
                    hits.append((a, e))
        bound = [(a, e) for a, e in hits if e["status"] == "bound"]
        lines.append(f"\n{cls} motif {mnum}  [{desc}]")
        lines.append(f"  members with a PDB: {len({a for a,_ in hits})} | "
                     f"structures: {len({e['pdb'] for _,e in hits})} | "
                     f"BOUND: {len({e['pdb'] for _,e in bound})}")
        for a, e in sorted(bound, key=lambda x: x[1]["pdb"]):
            lines.append(f"    * {e['pdb']}  [{e['detail']}]  {a}  {e['title']}")
            hit_rows.append([cls, mnum, a, e["pdb"], e["status"], e["detail"], e["title"]])
    report = "\n".join(lines)
    print("\n" + report, flush=True)
    (OUT / "report.txt").write_text(report + "\n")
    with open(OUT / "bound_structure_hits.csv", "w") as fh:
        w = csv.writer(fh)
        w.writerow(["class", "motif", "accession", "pdb", "status", "detail", "title"])
        w.writerows(hit_rows)
    print(f"\n[written] {OUT}/report.txt + bound_structure_hits.csv", flush=True)


if __name__ == "__main__":
    main()
