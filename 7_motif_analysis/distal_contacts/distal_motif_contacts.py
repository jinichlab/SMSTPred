#!/usr/bin/env python3
"""What do the distal motifs touch?  For every PDB entry of every carrier of a distal motif, locate the motif by its matched
sequence on the carrier chain, then list its heavy-atom contacts (<= CUT A) with (a) other polymer chains (protein of another
entity / copy of the carrier / nucleic acid), (b) small-molecule ligands and ions, and the surface it buries (SASA of the motif
residues in the carrier chain alone minus in the whole assembly).

  python distal_motif_contacts.py --dir <release>/motif_analysis --motifs DNA-4 RNA-5 ... --out OUT [--cif-cache DIR] [--max-per-motif 150]
Downloads missing mmCIF files from files.rcsb.org into the cache (public PDB entries).
"""
import argparse, json, os, re, sys, urllib.request, warnings
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path
import pandas as pd
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import three_to_index, index_to_one
warnings.filterwarnings("ignore")
CUT = 4.5
SAM_LIKE = {"SAM", "SAH", "SFG", "MTA", "DZF", "SX0", "EEM", "5AD", "SAE", "SMM", "AAT", "A5A", "SAI", "AMO"}
BUFFER = {"HOH", "GOL", "EDO", "PEG", "PGE", "PG4", "SO4", "PO4", "ACT", "FMT", "DMS", "TRS", "MES", "EPE", "CL", "NA", "BR", "IOD", "MPD", "BME", "PE4", "1PE", "P6G", "CIT", "TLA", "NO3", "IMD", "ACY", "MLI", "CAC"}
IONS = {"MG", "ZN", "MN", "CA", "K", "FE", "FE2", "CO", "NI", "CU", "CD", "SF4", "FES"}
AA = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","MSE"}
T2O = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M"}

def fetch(pdb, cache):
    p = cache / f"{pdb}.cif"
    if p.exists() and p.stat().st_size > 1000:
        return p
    try:
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb}.cif", p)
        return p
    except Exception as e:
        return None

def chain_info(meta):
    ch = {}
    for e in meta.get("polymer_entities") or []:
        t = e["entity_poly"]["rcsb_entity_polymer_type"]
        ids = e["rcsb_polymer_entity_container_identifiers"]
        refs = ids.get("reference_sequence_identifiers") or []
        acc = next((r["database_accession"] for r in refs if r["database_name"] == "UniProt"), "")
        for c in ids["auth_asym_ids"]:
            ch[c] = (t, acc, tuple(ids["auth_asym_ids"]))
    return ch

def work(args):
    row, cache, meta = args
    pdb = row["pdb"]; acc = row["accession"]; seq = row["matched_sequence"]
    f = fetch(pdb, cache)
    if f is None:
        return dict(motif=row["motif"], pdb=pdb, accession=acc, status="download_failed")
    try:
        s = MMCIFParser(QUIET=True).get_structure(pdb, str(f))
    except Exception as e:
        return dict(motif=row["motif"], pdb=pdb, accession=acc, status="parse_failed")
    model = s[0]
    cinfo = chain_info(meta)
    carriers = [c for c, v in cinfo.items() if v[1] == acc and c in model]
    out = dict(motif=row["motif"], pdb=pdb, accession=acc, carrier_in_class=row["in_class"], status="motif_not_resolved")
    for cid in carriers:
        chain = model[cid]
        res = [r for r in chain if r.id[0] == " " or r.get_resname() == "MSE"]
        res = [r for r in res if r.get_resname() in AA]
        sq = "".join(T2O[r.get_resname()] for r in res)
        i = sq.find(seq)
        if i < 0:
            continue
        mres = res[i:i + len(seq)]
        mset = set(id(r) for r in mres)
        atoms = [a for a in model.get_atoms()]
        ns = NeighborSearch(atoms)
        cont = defaultdict(set)   # category -> set of partner labels
        for r in mres:
            for a in r:
                for b in ns.search(a.coord, CUT):
                    rb = b.get_parent()
                    if rb is r:
                        continue
                    cb = rb.get_parent().id
                    rn = rb.get_resname().strip()
                    if cb == cid and rb.id[0] == " ":
                        continue  # intra-chain polymer contact
                    if rb.id[0] == "W" or rn == "HOH":
                        continue
                    if rb.id[0] == " " or rn in AA or rn in {"DA","DC","DG","DT","A","C","G","U","DU","N","DN"}:
                        t, acc2, ents = cinfo.get(cb, ("?", "", ()))
                        if t == "Protein":
                            lab = "protein:self-entity" if cb in cinfo[cid][2] else "protein:" + (acc2 or "other")
                        elif t in ("DNA", "RNA"):
                            lab = "nucleic-acid:" + t
                        else:
                            lab = t
                        cont[lab].add(f"{cb}:{rn}{rb.id[1]}")
                    else:
                        if rn in SAM_LIKE:
                            cont["ligand:SAM-like"].add(f"{cb}:{rn}")
                        elif rn in IONS:
                            cont["ion:" + rn].add(f"{cb}:{rn}")
                        elif rn in BUFFER:
                            cont["buffer"].add(f"{cb}:{rn}")
                        else:
                            cont["ligand:" + rn].add(f"{cb}:{rn}{rb.id[1]}")
        # buried surface by partners
        sr = ShrakeRupley()
        sr.compute(model, level="R")
        full = sum(r.sasa for r in mres)
        import copy
        alone = None
        try:
            m2 = copy.deepcopy(model)
            for c in list(m2.child_list):
                if c.id != cid:
                    m2.detach_child(c.id)
            for r in list(m2[cid]):
                if r.id[0] != " ":
                    m2[cid].detach_child(r.id)
            sr.compute(m2, level="R")
            mres2 = [m2[cid][r.id] for r in mres if r.id in m2[cid]]
            alone = sum(r.sasa for r in mres2)
        except Exception:
            pass
        out.update(status="ok", chain=cid, motif_range=f"{mres[0].id[1]}-{mres[-1].id[1]}", n_res=len(mres),
                   sasa_complex=round(full, 1), sasa_alone=None if alone is None else round(alone, 1),
                   buried_by_partners=None if alone is None else round(alone - full, 1),
                   contacts=json.dumps({k: sorted(v)[:6] for k, v in cont.items()}),
                   contact_classes=";".join(sorted(cont)))
        return out
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True); ap.add_argument("--motifs", nargs="+", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--cif-cache", default=None)
    ap.add_argument("--max-per-motif", type=int, default=150); ap.add_argument("--procs", type=int, default=24)
    a = ap.parse_args()
    d = Path(a.dir) / "structure_search"; out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    cache = Path(a.cif_cache) if a.cif_cache else out / "cif"; cache.mkdir(parents=True, exist_ok=True)
    meta = json.load(open(d / "pdb_entries_v2.json"))
    c = pd.read_csv(d / "motif_structure_candidates.csv")
    jobs = []
    for m in a.motifs:
        x = c[c.motif == m].copy()
        x["res"] = pd.to_numeric(x.resolution, errors="coerce")
        x = x.sort_values(["in_class", "res"], ascending=[False, True]).drop_duplicates("pdb").head(a.max_per_motif)
        for _, r in x.iterrows():
            jobs.append((r.to_dict(), cache, meta[r["pdb"]]))
    print(len(jobs), "jobs", flush=True)
    with Pool(a.procs) as p:
        res = p.map(work, jobs, chunksize=4)
    df = pd.DataFrame(res); df.to_csv(out / "contacts_per_structure.csv", index=False)
    print(df.groupby("motif").status.value_counts())
    ok = df[df.status == "ok"]
    lines = []
    for m, g in ok.groupby("motif"):
        n = len(g); cls = Counter()
        for s in g.contact_classes.fillna(""):
            for k in set(s.split(";")) - {""}:
                cls[k] += 1
        lines.append(f"{m}: {n} structures; buried by partners median {g.buried_by_partners.median()} A2 (max {g.buried_by_partners.max()}); contact classes: " + ", ".join(f"{k} {v}" for k, v in cls.most_common(14)))
    (out / "summary.txt").write_text("\n".join(lines) + "\n"); print("\n".join(lines))
main()
