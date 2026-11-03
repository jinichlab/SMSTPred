#!/usr/bin/env python3
"""
Stage 0 -- download Rhea's own classification of every reaction participant, for the substrate typing (substrate_binning.py).

  python fetch_rhea_compounds.py [--out rhea_compounds.tsv.gz]

One SPARQL query to https://sparql.rhea-db.org/sparql returns every compound Rhea knows (about 15,000) with its name and its family:

  Polypeptide      rh:GenericPolypeptide   a protein-bound residue, e.g. "L-histidyl-[protein]" (RHEA-COMP ids; 1,268 on 2026-11-02)
  Polynucleotide   rh:GenericPolynucleotide RNA / DNA-bound residues, e.g. "adenosine(37) in tRNA(Phe)" (712; every name carries an RNA or DNA word)
  Polymer          rh:Polymer              sugar chains, polyphosphate, teichoic acid ... (264)
  SmallMolecule    rh:SmallMolecule        ChEBI compounds (about 13,000)

The UniProt reaction text uses Rhea's compound names verbatim, so a participant is typed by an exact lookup of its name in this table: the family
comes from Rhea's own ontology, not from matching words inside a bracket. A polymer Rhea adds later appears with the next download. Names that are
not in the table (reactions UniProt annotated without a Rhea id) fall back to the word rules and are listed in the binning report.

Output: rhea_compounds.tsv.gz (accession, name, family) and rhea_compounds.json (query date, counts, sha256). Network: one request, a few seconds.
"""
import argparse
import gzip
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

ENDPOINT = "https://sparql.rhea-db.org/sparql"
QUERY = """
PREFIX rh: <http://rdf.rhea-db.org/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
SELECT DISTINCT ?acc ?name ?fam WHERE {
  ?c rh:accession ?acc ; rh:name ?name .
  FILTER(!STRSTARTS(?acc, "RHEA"))
  OPTIONAL { ?c rdfs:subClassOf ?fam . FILTER(?fam IN (rh:GenericPolypeptide, rh:GenericPolynucleotide, rh:Polymer, rh:SmallMolecule)) }
}"""
SHORT = {"GenericPolypeptide": "Polypeptide", "GenericPolynucleotide": "Polynucleotide", "Polymer": "Polymer", "SmallMolecule": "SmallMolecule"}


def fetch(retries=4):
    for attempt in range(retries):
        try:
            r = requests.get(ENDPOINT, params={"query": QUERY}, headers={"Accept": "text/tab-separated-values"}, timeout=180)
            r.raise_for_status()
            return r.text
        except requests.RequestException as e:
            wait = 5 * (attempt + 1)
            print(f"  attempt {attempt + 1} failed ({e}); retrying in {wait} s", file=sys.stderr)
            time.sleep(wait)
    sys.exit("the Rhea SPARQL endpoint could not be reached")


def unquote(v):
    v = v.strip()
    if v.startswith('"'):
        v = v[1:v.rindex('"')]
    return v.replace('\\"', '"').replace("\\\\", "\\")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "rhea_compounds.tsv.gz"))
    args = ap.parse_args()
    text = fetch()
    rows = {}
    for line in text.splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        acc, name, fam = unquote(parts[0]), unquote(parts[1]), parts[2].strip()
        fam = SHORT.get(fam.rsplit("/", 1)[-1], "Other") if fam else "Other"
        rows.setdefault((acc, name), fam)
    if len(rows) < 10000:
        sys.exit(f"only {len(rows)} compounds returned: the query result looks incomplete, nothing written")
    out = Path(args.out)
    with gzip.open(out, "wt") as f:
        f.write("accession\tname\tfamily\n")
        for (acc, name), fam in sorted(rows.items()):
            f.write(f"{acc}\t{name}\t{fam}\n")
    fam_counts = {}
    for fam in rows.values():
        fam_counts[fam] = fam_counts.get(fam, 0) + 1
    meta = {"endpoint": ENDPOINT, "downloaded": datetime.now(timezone.utc).isoformat(timespec="seconds"), "compounds": len(rows), "families": fam_counts,
            "sha256": hashlib.sha256(out.read_bytes()).hexdigest()}
    out.with_suffix("").with_suffix(".json").write_text(json.dumps(meta, indent=2))
    print(f"wrote {out}: {len(rows)} compounds {fam_counts}")


if __name__ == "__main__":
    main()
