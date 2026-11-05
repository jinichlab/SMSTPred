#!/usr/bin/env python3
"""Which proteins that carry a motif have an experimental structure with SAM (or an analogue) AND a substrate bound?

  python motif_structure_search.py --dir <release>/motif_analysis [--min-fold 2]

Step 1 of picking the structure that represents each motif, automatically (replaces the hand-picked panels of Figs 9/10).

Reads    <dir>/motif_specificity_matrix.csv   the motifs and their fold (fimo_cross_manifest.py); only motifs with ratio_vs_bg >= --min-fold are searched
         <dir>/fimo_cross/<class>.tsv          every FIMO hit (p < 1e-4) of that class's motifs on every protein of the corpus (run_fimo_cross.sh)
         <dir>/corpus.csv                      each protein's class
Queries  rest.uniprot.org   one streamed search: every UniProt entry with InterPro IPR029063 and a PDB cross-reference (its PDB ids)
         data.rcsb.org      GraphQL, 100 entries per request: method, resolution, the polymer entities (type, length, UniProt accession, chains) and the
                            small-molecule ligands with the chains they sit on.   (No structure file is downloaded; answers are cached in <dir>/structure_search/.)
Writes   <dir>/structure_search/
           pdb_entries_v2.json          cache of the RCSB answers (a re-run asks only for new entries)
           motif_structure_candidates.csv   one row per (motif, carrier, PDB entry): carrier class / reviewed / FIMO p-value / hit position, the entry's method and
                                            resolution, SAM-like ligands, substrate type and ligands, and `tier`
           motif_structure_summary.csv/.txt   per motif: carriers, carriers with a PDB, with SAM, with SAM + substrate (any carrier / carriers of the motif's own class)

Tiers of an entry for a carrier (the highest that holds):
  sam_substrate_on_chain   a SAM-like ligand AND a substrate, both on the carrier's own chain(s) (nucleic-acid / peptide substrates: present in the entry)
  sam_substrate            both are in the entry, but the SAM-like ligand is not on the carrier's chain (another protein of the complex holds it)
  sam_only                 SAM-like ligand, no substrate
  substrate_only           a substrate, no SAM-like ligand
  apo                      neither
SAM-like = SAM, SAH, SFG (sinefungin) and the analogues in SAM_LIKE. Substrate = a nucleic-acid polymer (DNA/RNA MTases), a peptide of at most PEPTIDE_MAX residues
that is not the carrier (protein MTases), or a small-molecule ligand of at least 100 Da that is not a cofactor of another enzyme, iron-sulfur cluster, unknown ligand, buffer or cryo-protectant (IGNORE; metabolite MTases).
This is an entry-level reading of the ligand lists: it does not yet check that the motif lies in the resolved part of the structure or near the ligands (that needs the coordinates).
"""
import argparse
import json
import re
import time
from pathlib import Path

import pandas as pd
import requests

CLASSES = ["DNA", "RNA", "Protein", "Metabolite"]
SAM_LIKE = {"SAM", "SAH", "SFG", "MTA", "DZF", "SX0", "EEM", "5AD", "SAE", "SMM", "AAT", "A5A", "SAI", "AMO", "ADE_SAM"}
IGNORE = {"HOH", "GOL", "EDO", "PEG", "PG4", "PG0", "1PE", "P6G", "SO4", "PO4", "CL", "NA", "MG", "K", "ZN", "CA", "MN", "NI", "FE", "FE2", "CO", "CD", "CU",
          "ACT", "DMS", "IOD", "BR", "FMT", "MPD", "TRS", "EPE", "BME", "ACE", "NO3", "CAC", "IMD", "BCT", "MES", "FLC", "GG", "PGE", "PEU", "OXY", "SCN",
          "DTT", "TLA", "CIT", "MLI", "NH4", "AZI", "HEZ", "1PG", "2PE", "PGO", "PG6", "BU3", "15P", "PE4", "PE5", "DIO", "EOH", "IPA", "MOH", "TAM", "BTB", "ACY"}
# not a substrate either: iron-sulfur clusters, unknown atoms / ligands, nucleotide and redox cofactors of other enzymes (a ligand under MIN_LIGAND_DA is an ion or a solvent)
IGNORE |= {"FES", "SF4", "F3S", "3FS", "UNX", "UNL", "ATP", "ADP", "AMP", "ANP", "GTP", "GDP", "GNP", "NAD", "NAP", "NDP", "NAI", "FAD", "FMN", "HEM", "HEC", "PLP", "COA", "ACO",
           "PNS", "TPP", "UDP", "UTP", "U5P", "CTP", "CDP"}
MIN_LIGAND_DA = 100
PEPTIDE_MAX = 30
GRAPHQL = "https://data.rcsb.org/graphql"
QUERY = """query($ids:[String!]!){ entries(entry_ids:$ids){ rcsb_id exptl{method} rcsb_entry_info{resolution_combined} struct{title}
 polymer_entities{ entity_poly{rcsb_entity_polymer_type rcsb_sample_sequence_length}
   rcsb_polymer_entity_container_identifiers{auth_asym_ids reference_sequence_identifiers{database_name database_accession}} }
 nonpolymer_entities{ nonpolymer_comp{chem_comp{id name type formula_weight}} rcsb_nonpolymer_entity_container_identifiers{auth_asym_ids} } } }"""


def get(url, **kw):
    """GET with retries; a rate limit (429) or server error waits as long as the server asks (Retry-After) or an increasing time, up to about 10 minutes in total."""
    for attempt in range(8):
        try:
            r = requests.get(url, timeout=120, **kw)
            if r.status_code == 429 or r.status_code >= 500:
                wait = int(r.headers.get("Retry-After", 0) or 0) or 30 * (attempt + 1)
                if attempt == 7:
                    r.raise_for_status()
                time.sleep(min(wait, 120))
                continue
            r.raise_for_status()
            return r
        except requests.exceptions.RequestException:
            if attempt == 7:
                raise
            time.sleep(10 * (attempt + 1))


def uniprot_pdb_map():
    """({UniProt accession: [PDB ids]}, {accession: reviewed / protein name / organism}) for every IPR029063 entry that has a PDB cross-reference.
    Read in pages of 500 by following the cursor link (a single streamed request of this size times out on busy days)."""
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {"query": "(xref:interpro-IPR029063) AND (database:pdb)", "fields": "accession,reviewed,protein_name,organism_name,xref_pdb", "format": "tsv", "size": 500}
    out, names, page = {}, {}, 0
    while url:
        r = get(url, params=params)
        params = None                       # the next link carries the query and the cursor
        lines = r.text.splitlines()
        for line in lines[1:]:
            acc, reviewed, name, org, pdbs = (line.split("\t") + [""] * 5)[:5]
            ids = [x for x in pdbs.split(";") if x]
            if ids:
                out[acc] = ids
                names[acc] = {"reviewed": reviewed == "reviewed", "protein_name": name, "organism": org}
        page += 1
        if page % 20 == 0:
            print(f"[uniprot] {page} pages, {len(out)} entries with PDB", flush=True)
        url = (r.links.get("next") or {}).get("url")
    return out, names


def rcsb_entries(ids, cache):
    """Fill `cache` ({PDB id: RCSB record}) for the ids it lacks, 100 per GraphQL request."""
    todo = [i for i in ids if i not in cache]
    for k in range(0, len(todo), 100):
        batch = todo[k:k + 100]
        for attempt in range(5):
            try:
                r = requests.post(GRAPHQL, json={"query": QUERY, "variables": {"ids": batch}}, timeout=180)
                r.raise_for_status()
                data = r.json()["data"]["entries"]
                break
            except Exception:
                if attempt == 4:
                    raise
                time.sleep(3 * (attempt + 1))
        got = {e["rcsb_id"]: e for e in data if e}
        for i in batch:
            cache[i] = got.get(i)       # None: RCSB does not know the entry (obsolete); kept so it is not asked again
        print(f"[rcsb] {min(k + 100, len(todo))}/{len(todo)}", flush=True)
    return cache


def classify(entry, accession):
    """The reading of one RCSB record for one carrier. Returns a dict: sam, substrate_type, substrate, tier, method, resolution, title."""
    chains, nucleic, peptides = set(), False, []
    for p in entry.get("polymer_entities") or []:
        ids = p["rcsb_polymer_entity_container_identifiers"]
        refs = {x["database_accession"] for x in (ids.get("reference_sequence_identifiers") or []) if x.get("database_name") == "UniProt"}
        typ, ln = p["entity_poly"]["rcsb_entity_polymer_type"], p["entity_poly"]["rcsb_sample_sequence_length"] or 0
        if accession in refs:
            chains |= set(ids.get("auth_asym_ids") or [])
        elif typ in ("DNA", "RNA", "NA-hybrid"):
            nucleic = True
        elif typ == "Protein" and ln <= PEPTIDE_MAX:
            peptides.append(ln)
    sam, sam_on_chain, small, small_on_chain = set(), False, set(), False
    for n in entry.get("nonpolymer_entities") or []:
        cc = n["nonpolymer_comp"]["chem_comp"]
        comp = cc["id"]
        on = bool(chains & set(n["rcsb_nonpolymer_entity_container_identifiers"].get("auth_asym_ids") or []))
        if comp in SAM_LIKE:
            sam.add(comp)
            sam_on_chain |= on
        elif comp not in IGNORE and (cc.get("formula_weight") or 999) >= MIN_LIGAND_DA:
            small.add(comp)
            small_on_chain |= on
    sub_type = [t for t, present in (("nucleic_acid", nucleic), ("peptide", bool(peptides)), ("small_molecule", bool(small))) if present]
    substrate = bool(sub_type)
    if sam and substrate:
        tier = "sam_substrate_on_chain" if (sam_on_chain and (nucleic or peptides or small_on_chain)) else "sam_substrate"
    elif sam:
        tier = "sam_only"
    elif substrate:
        tier = "substrate_only"
    else:
        tier = "apo"
    res = (entry.get("rcsb_entry_info") or {}).get("resolution_combined") or [None]
    return {"sam": "+".join(sorted(sam)), "substrate_type": "+".join(sub_type), "substrate": "+".join(sorted(small)), "tier": tier,
            "method": ((entry.get("exptl") or [{}])[0]).get("method", ""), "resolution": res[0], "title": (entry.get("struct") or {}).get("title", "")[:100]}


TIERS = ["sam_substrate_on_chain", "sam_substrate", "sam_only", "substrate_only", "apo"]


def carriers(d, kept):
    """FIMO hits of the kept motifs: best (lowest p) hit per (motif, protein)."""
    frames = []
    for cls in CLASSES:
        f = d / "fimo_cross" / f"{cls}.tsv"
        if not f.exists():
            raise SystemExit(f"{f} not found: run 7_motif_analysis/run_fimo_cross.sh first")
        h = pd.read_csv(f, sep="\t", comment="#").dropna(subset=["sequence_name"])
        h["motif"] = cls + "-" + h["motif_id"].astype(str).str.split("-").str[0]
        frames.append(h[h["motif"].isin(kept)])
    h = pd.concat(frames, ignore_index=True)
    h = h.sort_values("p-value").drop_duplicates(["motif", "sequence_name"])
    return h.rename(columns={"sequence_name": "accession", "p-value": "fimo_p"})[["motif", "accession", "start", "stop", "fimo_p", "matched_sequence"]]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="<release>/motif_analysis")
    ap.add_argument("--min-fold", type=float, default=2.0)
    a = ap.parse_args()
    d = Path(a.dir)
    out = d / "structure_search"
    out.mkdir(exist_ok=True)
    spec = pd.read_csv(d / "motif_specificity_matrix.csv")
    kept = spec[spec["ratio_vs_bg"] >= a.min_fold]
    print(f"{len(kept)} of {len(spec)} motifs reach {a.min_fold:g}x", flush=True)
    hits = carriers(d, set(kept["motif"]))
    corpus = pd.read_csv(d / "corpus.csv", usecols=["Entry", "class"], keep_default_na=False, low_memory=False)
    hits = hits.merge(corpus.rename(columns={"Entry": "accession", "class": "carrier_class"}), on="accession", how="left")
    hits["in_class"] = hits["motif"].str.split("-").str[0] == hits["carrier_class"]
    print(f"{len(hits):,} (motif, protein) carriers over {hits.accession.nunique():,} proteins", flush=True)

    acc2pdb, names = uniprot_pdb_map()
    hits = hits[hits["accession"].isin(acc2pdb)]
    pairs = [(m, acc, p) for m, acc in zip(hits["motif"], hits["accession"]) for p in acc2pdb[acc]]
    pdbs = sorted({p for _, _, p in pairs})
    print(f"{hits.accession.nunique():,} carriers have a PDB cross-reference: {len(pdbs):,} entries", flush=True)

    cache_file = out / "pdb_entries_v2.json"       # v2: records carry formula_weight
    cache = json.load(open(cache_file)) if cache_file.exists() else {}
    rcsb_entries(pdbs, cache)
    json.dump(cache, open(cache_file, "w"))

    rows = []
    byhit = hits.set_index(["motif", "accession"])
    for m, acc, p in pairs:
        e = cache.get(p)
        if not e:
            continue
        h = byhit.loc[(m, acc)]
        rows.append({"motif": m, "accession": acc, "carrier_class": h["carrier_class"], "in_class": bool(h["in_class"]), "reviewed": names[acc]["reviewed"], "protein_name": names[acc]["protein_name"], "organism": names[acc]["organism"],
                     "fimo_p": h["fimo_p"], "start": h["start"], "stop": h["stop"], "matched_sequence": h["matched_sequence"], "pdb": p, **classify(e, acc)})
    cand = pd.DataFrame(rows)
    cand["tier_rank"] = cand["tier"].map({t: i for i, t in enumerate(TIERS)})
    cand = cand.sort_values(["motif", "tier_rank", "in_class", "resolution", "fimo_p"], ascending=[True, True, False, True, True]).drop(columns="tier_rank")
    cand.to_csv(out / "motif_structure_candidates.csv", index=False)

    n_carr = hits.groupby("motif").accession.nunique()
    allc = carriers(d, set(kept["motif"])).groupby("motif").accession.nunique()
    srows = []
    for m in kept["motif"]:
        c = cand[cand["motif"] == m]
        both = c[c["tier"].isin(["sam_substrate_on_chain", "sam_substrate"])]
        srows.append({"motif": m, "consensus": kept.set_index("motif").loc[m, "consensus"], "fold": kept.set_index("motif").loc[m, "ratio_vs_bg"],
                      "carriers": int(allc.get(m, 0)), "with_pdb": int(c.accession.nunique()), "with_sam": int(c[c.sam != ""].accession.nunique()),
                      "with_sam_and_substrate": int(both.accession.nunique()), "with_sam_and_substrate_in_class": int(both[both.in_class].accession.nunique()),
                      "sam_substrate_on_chain": int(c[c.tier == "sam_substrate_on_chain"].accession.nunique()),
                      "best": (both.iloc[0]["accession"] + " " + both.iloc[0]["pdb"]) if len(both) else ""})
    s = pd.DataFrame(srows)
    s.to_csv(out / "motif_structure_summary.csv", index=False)
    (out / "motif_structure_summary.txt").write_text(s.to_string(index=False) + "\n")
    print(s.to_string(index=False))
    print(f"\nwrote {out}/motif_structure_candidates.csv, motif_structure_summary.csv/.txt")


if __name__ == "__main__":
    main()
