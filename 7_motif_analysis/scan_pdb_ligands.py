#!/usr/bin/env python3
"""Characterize every PDB in the motif carrier table for SAM cofactor + substrate.

Reads the distinct PDB IDs from motif_structure_carriers.csv and queries RCSB
(GraphQL data API; login node only) for each structure's polymer content
(DNA / RNA / protein chain counts) and its bound non-polymer ligands. Classifies:

  cofactor : a SAM-family cofactor is bound (SAM / SAH / sinefungin / SAM-analog)
  substrate: class-dependent, applied when joined to a motif --
     DNA motif  -> a DNA chain present
     RNA motif  -> an RNA chain present
     Protein    -> a second protein chain (candidate substrate/partner)
     Metabolite      -> a non-cofactor, non-additive small molecule present

Outputs:
  pdb_ligand_characterization.csv   one row per PDB: polymer counts, cofactor, ligands
  motif_structure_ligands.csv       carrier table + per-row has_cofactor/has_substrate
  scan_pdb_ligands.txt              per-motif: how many carriers give SAM + substrate
"""
import json
import time
import urllib.request
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
GRAPHQL = "https://data.rcsb.org/graphql"
BATCH = 50

# SAM-family cofactor ligand codes (verified SAM/SAH/analogs); we ALSO flag any
# ligand whose name contains "adenosyl", to catch analogs not in this list.
COFACTOR = {"SAM", "SAH", "SFG", "0XU", "NEA", "QVR", "MTA", "EEM", "SX0"}
# crystallisation additives / buffer / ions -- never a substrate
ADDITIVE = {"HOH", "GOL", "EDO", "ACT", "PO4", "SO4", "NA", "CL", "K", "MG", "ZN",
            "CA", "MN", "NI", "FE", "CO", "CD", "HG", "BME", "DMS", "IOD", "BR",
            "ACE", "FMT", "ADP", "ATP", "GDP", "GTP", "EPE", "HF2", "UNX", "PEG",
            "PG4", "PGE", "1PE", "MPD", "TRS", "IMD", "CIT", "FLC", "MES", "BTB",
            "TLA", "MLI", "SCN", "NO3", "NH4", "PGO", "BCT", "CAC", "DTT", "TCE"}


def query(ids):
    q = ('{ entries(entry_ids: %s) { rcsb_id '
         'rcsb_entry_info { polymer_entity_count_DNA polymer_entity_count_RNA '
         'polymer_entity_count_protein selected_polymer_entity_types } '
         'nonpolymer_entities { pdbx_entity_nonpoly { comp_id name } } } }'
         % json.dumps(ids))
    for attempt in range(4):
        try:
            req = urllib.request.Request(GRAPHQL, data=json.dumps({"query": q}).encode(),
                                         headers={"Content-Type": "application/json",
                                                  "User-Agent": "mtase-ligand-scan"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)["data"]["entries"] or []
        except Exception as e:
            if attempt == 3:
                print(f"  [warn] batch failed: {e}"); return []
            time.sleep(3 * (attempt + 1))


def main():
    carriers = pd.read_csv(HERE / "motif_structure_carriers.csv")
    pdbs = sorted({p for s in carriers.pdbs.dropna() for p in s.split(";")})
    print(f"[scan] {len(pdbs)} distinct PDB structures")

    recs = {}
    for i in range(0, len(pdbs), BATCH):
        for e in query([p.upper() for p in pdbs[i:i + BATCH]]):
            ei = e["rcsb_entry_info"]
            ligs = [(n["pdbx_entity_nonpoly"]["comp_id"],
                     (n["pdbx_entity_nonpoly"].get("name") or "").lower())
                    for n in (e["nonpolymer_entities"] or [])]
            cof = sorted({c for c, nm in ligs if c in COFACTOR or "adenosyl" in nm})
            others = sorted({c for c, nm in ligs
                             if c not in COFACTOR and "adenosyl" not in nm and c not in ADDITIVE})
            recs[e["rcsb_id"].upper()] = {
                "pdb": e["rcsb_id"].upper(),
                "poly_types": ei["selected_polymer_entity_types"],
                "n_dna": ei["polymer_entity_count_DNA"],
                "n_rna": ei["polymer_entity_count_RNA"],
                "n_protein": ei["polymer_entity_count_protein"],
                "cofactor": ";".join(cof), "has_cofactor": bool(cof),
                "small_mol": ";".join(others), "n_small_mol": len(others),
            }
        print(f"  {min(i+BATCH, len(pdbs))}/{len(pdbs)}", flush=True)

    lig = pd.DataFrame(recs.values())
    lig.to_csv(HERE / "pdb_ligand_characterization.csv", index=False)
    miss = set(p.upper() for p in pdbs) - set(recs)
    print(f"[scan] characterized {len(lig)}; {len(miss)} unresolved (obsolete/superseded)")
    print(f"  with SAM-family cofactor : {lig.has_cofactor.sum()}")
    print(f"  with DNA chain           : {(lig.n_dna>0).sum()}")
    print(f"  with RNA chain           : {(lig.n_rna>0).sum()}")
    print(f"  with >=2 protein chains  : {(lig.n_protein>=2).sum()}")
    print(f"  with a small-molecule    : {(lig.n_small_mol>0).sum()}")

    # ---- join to carriers; substrate is class-dependent ----
    L = lig.set_index("pdb")
    def characterize(row):
        # use the carrier's FIRST listed PDB as the representative structure
        pdb = str(row.pdbs).split(";")[0].upper()
        if pdb not in L.index:
            return pd.Series({"rep_pdb": pdb, "has_cofactor": False,
                              "has_substrate": False, "substrate_kind": "?", "poly_types": "?"})
        r = L.loc[pdb]; cls = row.motif_class
        if cls == "DNA":
            sub, kind = r.n_dna > 0, "DNA"
        elif cls == "RNA":
            sub, kind = r.n_rna > 0, "RNA"
        elif cls == "Protein":
            sub, kind = r.n_protein >= 2, "partner protein"
        else:
            sub, kind = r.n_small_mol > 0, r.small_mol
        return pd.Series({"rep_pdb": pdb, "has_cofactor": bool(r.has_cofactor),
                          "has_substrate": bool(sub), "substrate_kind": kind if sub else "",
                          "poly_types": r.poly_types})
    out = pd.concat([carriers, carriers.apply(characterize, axis=1)], axis=1)
    out["cofactor_and_substrate"] = out.has_cofactor & out.has_substrate
    out.to_csv(HERE / "motif_structure_ligands.csv", index=False)

    # ---- per-motif availability (in-class carriers with BOTH) ----
    rep = ["Structure availability per specific motif "
           "(in-class carriers with a bound SAM cofactor AND substrate)\n"]
    rep.append(f"{'motif':10} {'consensus':18} {'in-class w/':11} {'best carrier (SAM + substrate)':52}")
    rep.append(f"{'':10} {'':18} {'SAM+subst':11}")
    for motif, g in out.groupby("motif"):
        cons = g.consensus.iloc[0]
        good = g[g.in_class & g.cofactor_and_substrate].sort_values(
            ["reviewed", "p_value"], ascending=[False, True])
        n = len(good)
        if n:
            b = good.iloc[0]
            best = (f"{b.accession} {str(b.protein_name)[:24]:24} {b.rep_pdb} "
                    f"[cof {L.loc[b.rep_pdb,'cofactor']}; sub {b.substrate_kind}] @{b.fimo_start}-{b.fimo_stop}")
        else:
            # fall back: any in-class with cofactor only, or note gap
            cofonly = g[g.in_class & g.has_cofactor]
            best = (f"none with both; {len(cofonly)} in-class have SAM only"
                    if len(cofonly) else "!! none in-class with SAM at all")
        rep.append(f"{motif:10} {cons:18} {n:>11} {best}")
    txt = "\n".join(rep)
    (HERE / "scan_pdb_ligands.txt").write_text(txt + "\n")
    print("\n" + txt)
    print("\n[wrote] pdb_ligand_characterization.csv + motif_structure_ligands.csv + scan_pdb_ligands.txt")


if __name__ == "__main__":
    main()
