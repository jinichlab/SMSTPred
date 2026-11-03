#!/usr/bin/env python3
"""Pick, for every motif, the experimental structure that represents it -- by rule, from coordinates -- and measure where the motif sits.

  python motif_structure_pick.py --dir <release>/motif_analysis [--max-try 40] [--max-pass 8] [--max-mb 60]

Step 2 of the automatic motif-structure system (step 1: motif_structure_search.py). For each motif, in this order, the carriers of motif_structure_candidates.csv that have SAM(-like)
and a substrate of the motif class's kind in the entry (DNA/RNA motifs: a nucleic acid; Protein: a peptide; Metabolite: a small-molecule ligand) are evaluated, in that order, until --max-pass of them have passed (or --max-try were tried), and then RANKED (rank_key): own-kind substrate, own-class carrier, the placement most of
the motif's structures show, reviewed, intact SAM/sinefungin, resolution, FIMO p-value. Rank 1 is drawn; all passing structures are in motif_pick_alternatives.csv. For each:
  (order of evaluation: own substrate kind, own class, reviewed, best FIMO p-value, best resolution)
  for each, the mmCIF file is downloaded (files.rcsb.org; kept in <dir>/structure_search/cif/; entries over --max-mb are skipped),
  the motif is LOCATED by searching the FIMO matched_sequence in the protein chains (author numbering; one mismatch allowed; a motif that is not in the resolved part of the structure fails),
  the SAM-like ligand and the substrate are identified (nucleic-acid chains / peptide chains of at most 30 residues / the small-molecule copy nearest the SAM-like ligand) and
  the minimum heavy-atom distances motif -> SAM-like ligand and motif -> substrate are measured, and substrate -> SAM-like ligand.
A structure passes when the SAM-like ligand and the substrate sit on this enzyme (motif -> each within GATE = 28 A: a ligand on another subunit of a big complex fails) and the substrate is
at the active site (substrate -> SAM-like ligand within 6 A for a small molecule, 8 A for nucleic acid / peptide).
Placement of the motif (as in the archived Table 2): "near substrate/SAM", "near SAM", "near substrate" (a distance of at most NEAR = 4 A: direct contact) or "distal" (neither).

Writes <dir>/structure_search/
  motif_placements.csv       one row per motif that has a passing structure: the columns of the archived structures/motif_placements.csv (motif, consensus, class, protein, organism, accession, pdb, method,
                             cofactor, substrate, chain, range, matched_seq, d_SAM, d_substrate, placement, folder) plus d_substrate_SAM, resolution, in_class, reviewed, tried, substrate_kind_fallback
  motif_pick_log.csv         every structure tried and why it failed or passed
  motif_pick_alternatives.csv   every passing structure of every motif with its rank (1 = the one drawn)
  structures/<Class>/<class>_motif<N>_<PDB>/   subset PDB (enzyme + SAM + substrate) with view.cxc (ChimeraX) and view.pml (PyMOL) scripts, styled as the archived panels
A motif whose carriers hold no entry with its class's substrate kind that passes falls back to the other kinds (substrate_kind_fallback = True). Motifs without a passing structure are listed in the log: they are not drawn.
"""
import argparse
import math
import sys
import time
from pathlib import Path

import gemmi
import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "structures"))
import make_context_views as mc  # noqa: E402   (subset(), cxc(), pml(): the styling of the archived panels)

CLASS_HEX = {"DNA": "6a51a3", "RNA": "cb181d", "Protein": "238b45", "Metabolite": "e6550d"}
EXPECT = {"DNA": "nucleic_acid", "RNA": "nucleic_acid", "Protein": "peptide", "Metabolite": "small_molecule"}
KIND = {"nucleic_acid": "nucleic", "peptide": "partner", "small_molecule": "ligand"}
COF = set(mc.COF_LIGS) | {"SAM", "SAH", "SFG", "MTA", "DZF", "SX0", "EEM", "5AD", "SAE", "SMM", "AAT", "A5A", "SAI", "AMO"}
NEAR = 4.0
GATE = 28.0
SUB_SAM_MAX = {"small_molecule": 6.0, "nucleic_acid": 8.0, "peptide": 8.0}
PEPTIDE_MAX = 30


def download(pdb, cif_dir, max_mb):
    """The mmCIF of an entry from files.rcsb.org into cif_dir; None if it is bigger than max_mb or unavailable."""
    f = cif_dir / f"{pdb}.cif"
    if f.exists():
        return f
    url = f"https://files.rcsb.org/download/{pdb}.cif"
    for attempt in range(4):
        try:
            with requests.get(url, stream=True, timeout=120) as r:
                r.raise_for_status()
                size = int(r.headers.get("Content-Length") or 0)
                if size > max_mb * 1e6:
                    return None
                tmp = f.with_suffix(".part")
                n = 0
                with open(tmp, "wb") as fh:
                    for chunk in r.iter_content(1 << 20):
                        n += len(chunk)
                        if n > max_mb * 1e6:
                            tmp.unlink()
                            return None
                        fh.write(chunk)
                tmp.replace(f)
                return f
        except Exception:
            time.sleep(3 * (attempt + 1))
    return None


def one(res):
    t = gemmi.find_tabulated_residue(res.name)
    return t.one_letter_code.upper() if t and t.is_amino_acid() else ""


def heavy(residues):
    return [(a.pos.x, a.pos.y, a.pos.z) for r in residues for a in r if a.element.name != "H"]


def mindist(A, B):
    return min((math.dist(p, q) for p in A for q in B), default=float("inf"))


def locate_all(model, matched):
    """Every protein chain that holds the motif: [(chain name, first residue number, last residue number, residues)]. One mismatch is allowed (a construct differs from the carrier)."""
    found = []
    for chain in model:
        aa = [r for r in chain if one(r)]
        seq = "".join(one(r) for r in aa)
        i = seq.find(matched)
        if i < 0:
            for j in range(len(seq) - len(matched) + 1):
                if sum(x != y for x, y in zip(seq[j:j + len(matched)], matched)) <= 1:
                    i = j
                    break
        if i >= 0:
            res = aa[i:i + len(matched)]
            found.append((chain.name, res[0].seqid.num, res[-1].seqid.num, res))
    return found


def is_nucleic(chain):
    return any((t := gemmi.find_tabulated_residue(r.name)) and t.is_nucleic_acid() for r in chain)


def evaluate(cif, matched, kind, sub_codes, carrier_chains=None):
    """Locate the motif and measure it in one structure. Returns (result dict | None, reason)."""
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    model = st[0]
    locs = locate_all(model, matched)
    if not locs:
        return None, "motif not found in the resolved protein chains"
    cof_res = [(ch.name, r) for ch in model for r in ch if r.name in COF]
    if not cof_res:
        return None, "no SAM-like ligand in the coordinates"
    best = None
    for chain, lo, hi, mres in locs:
        matoms = heavy(mres)
        # the SAM-like ligand copy nearest the motif
        cof = min(cof_res, key=lambda cr: mindist(matoms, heavy([cr[1]])))
        cof_at = heavy([cof[1]])
        if kind == "nucleic_acid":
            subs = [ch for ch in model if is_nucleic(ch)]
            sub_at = heavy([r for ch in subs for r in ch])
            keep = {chain, cof[0]} | {ch.name for ch in subs if mindist(cof_at, heavy(list(ch))) <= 12}
            sub_chains = [ch.name for ch in subs if ch.name in keep]
            sub_lig = None
        elif kind == "peptide":
            pep = [ch for ch in model if ch.name != chain and 0 < sum(1 for r in ch if one(r)) <= PEPTIDE_MAX and not is_nucleic(ch)]
            sub_at = heavy([r for ch in pep for r in ch if one(r)])
            sub_chains = [ch.name for ch in pep if mindist(cof_at, heavy([r for r in ch if one(r)])) <= 12]
            keep = {chain, cof[0]} | set(sub_chains)
            sub_lig = None
        else:
            cand = [(ch.name, r) for ch in model for r in ch if r.name in sub_codes]
            if not cand:
                continue
            c = min(cand, key=lambda cr: mindist(cof_at, heavy([cr[1]])))
            sub_at = heavy([c[1]])
            keep = {chain, cof[0], c[0]}
            sub_chains, sub_lig = [], c[1].name
        if not sub_at:
            continue
        d_cof, d_sub, d_ss = mindist(matoms, cof_at), mindist(matoms, sub_at), mindist(cof_at, sub_at)
        r = dict(st=st, model=model, chain=chain, lo=lo, hi=hi, d_cof=d_cof, d_sub=d_sub, d_ss=d_ss, keep=keep, sub_chains=sub_chains, sub_lig=sub_lig,
                 kind=KIND[kind], cof_name=cof[1].name)
        if best is None or (d_cof + d_sub) < (best["d_cof"] + best["d_sub"]):
            best = r
    if best is None:
        return None, "no substrate of the expected kind in the coordinates"
    if best["d_cof"] > GATE or best["d_sub"] > GATE:
        return None, f"SAM {best['d_cof']:.0f} A / substrate {best['d_sub']:.0f} A from the motif: a ligand sits on another subunit"
    if best["d_ss"] > SUB_SAM_MAX[kind]:
        return None, f"substrate {best['d_ss']:.1f} A from SAM: not at the active site"
    return best, "ok"


def rank_key(p, majority):
    """Sort key (smaller = better) of one passing structure `p` (dict: kind_fallback, in_class, placement, reviewed, cofactor, resolution, fimo_p); `majority` = the placement most
    of the motif's passing structures show. Compared in this order, each criterion deciding only among ties of the ones above:
      1 a substrate of the motif class's own kind (not a fallback)   2 a carrier of the motif's own class   3 the placement most structures of the motif show (the typical one)
      4 a reviewed (Swiss-Prot) carrier                                5 an intact SAM or sinefungin rather than the product SAH     6 resolution (0.5 A steps)   7 the FIMO p-value"""
    res = p["resolution"] if p["resolution"] == p["resolution"] and p["resolution"] is not None else 99.0      # NaN / None: unknown, last
    return (bool(p["kind_fallback"]), not p["in_class"], p["placement"] != majority, not p["reviewed"], not any(c in str(p["cofactor"]) for c in ("SAM", "SFG")),
            round(float(res) * 2) / 2, float(p["fimo_p"]))


def majority_placement(passing):
    """The placement shown by most of the structures (own-class, expected-kind ones if there are any); ties: the one nearer the active site."""
    order = ["near substrate/SAM", "near SAM", "near substrate", "distal"]
    pool = [p for p in passing if p["in_class"] and not p["kind_fallback"]] or passing
    counts = pd.Series([p["placement"] for p in pool]).value_counts()
    return min(counts.index, key=lambda k: (-counts[k], order.index(k)))


def placement(d_cof, d_sub):
    near_c, near_s = d_cof <= NEAR, d_sub <= NEAR
    return "near substrate/SAM" if near_c and near_s else "near SAM" if near_c else "near substrate" if near_s else "distal"


def method_of(cif):
    v = gemmi.cif.read(str(cif)).sole_block().find_value("_exptl.method") or "?"
    return str(v).strip().strip("'\"").replace("X-RAY DIFFRACTION", "X-ray").replace("ELECTRON MICROSCOPY", "cryo-EM").replace("SOLUTION NMR", "NMR")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="<release>/motif_analysis")
    ap.add_argument("--max-try", type=int, default=40)
    ap.add_argument("--max-pass", type=int, default=8, help="stop looking when this many structures of a motif have passed")
    ap.add_argument("--max-mb", type=float, default=60.0)
    a = ap.parse_args()
    d = Path(a.dir)
    out = d / "structure_search"
    cif_dir = out / "cif"
    cif_dir.mkdir(exist_ok=True)
    cand = pd.read_csv(out / "motif_structure_candidates.csv")
    cand = cand[cand["tier"].isin(["sam_substrate_on_chain", "sam_substrate"])]
    spec = pd.read_csv(d / "motif_specificity_matrix.csv").set_index("motif")
    rows, log, alts = [], [], []
    for motif in sorted(cand["motif"].unique(), key=lambda m: (m.split("-")[0], int(m.split("-")[1]))):
        cls, num = motif.split("-")
        c = cand[cand["motif"] == motif].copy()
        c["expected_kind"] = c["substrate_type"].str.contains(EXPECT[cls])
        # carriers whose entry holds the substrate kind of the motif's class come first; a motif that finds no passing structure that way falls back to any substrate kind (flagged)
        c = c.sort_values(["expected_kind", "in_class", "reviewed", "fimo_p", "resolution"], ascending=[False, False, False, True, True])
        passing, tried = [], 0
        for _, r in c.iterrows():
            if tried >= a.max_try or len(passing) >= a.max_pass:
                break
            tried += 1
            kind = EXPECT[cls] if r["expected_kind"] else next(k for k in ("nucleic_acid", "peptide", "small_molecule") if k in r["substrate_type"])
            cif = download(r["pdb"], cif_dir, a.max_mb)
            if cif is None:
                log.append({"motif": motif, "accession": r["accession"], "pdb": r["pdb"], "result": f"not downloaded (over {a.max_mb:g} MB or unavailable)"})
                continue
            ev, why = evaluate(cif, str(r["matched_sequence"]), kind, set(str(r["substrate"]).split("+")) if isinstance(r["substrate"], str) else set())
            log.append({"motif": motif, "accession": r["accession"], "pdb": r["pdb"], "in_class": r["in_class"], "kind": kind, "result": why})
            if ev is not None:
                passing.append({"r": r, "cif": cif, "ev": ev, "kind": kind, "kind_fallback": kind != EXPECT[cls], "in_class": bool(r["in_class"]), "reviewed": bool(r["reviewed"]),
                                "placement": placement(ev["d_cof"], ev["d_sub"]), "cofactor": r["sam"], "resolution": r["resolution"], "fimo_p": r["fimo_p"]})
        if not passing:
            print(f"  [NONE] {motif}: no passing structure among {tried} tried ({len(c)} eligible)", flush=True)
            continue
        maj = majority_placement(passing)
        passing.sort(key=lambda q: rank_key(q, maj))
        for rank, q in enumerate(passing, 1):
            alts.append({"motif": motif, "rank": rank, "chosen": rank == 1, "accession": q["r"]["accession"], "pdb": q["r"]["pdb"], "placement": q["placement"], "majority_placement": maj,
                         "in_class": q["in_class"], "reviewed": q["reviewed"], "kind_fallback": q["kind_fallback"], "cofactor": q["cofactor"], "resolution": q["resolution"],
                         "fimo_p": q["fimo_p"], "d_SAM": round(q["ev"]["d_cof"], 1), "d_substrate": round(q["ev"]["d_sub"], 1)})
        best = passing[0]
        r, cif, ev, kind = best["r"], best["cif"], best["ev"], best["kind"]
        place = best["placement"]
        pdb = r["pdb"]
        folder = f"{cls.lower()}_motif{num}_{pdb}"
        fd = out / "structures" / cls / folder
        fd.mkdir(parents=True, exist_ok=True)
        ev["st"].remove_alternative_conformations()
        tmp = fd / f"_full_{pdb}.pdb"
        ev["st"].write_pdb(str(tmp))
        rng = f"{ev['lo']}-{ev['hi']}"
        sub_short = {"nucleic": cls if cls in ("DNA", "RNA") else "nucleic acid", "partner": "peptide"}.get(ev["kind"], ev["sub_lig"])
        cap = f"{motif} {r['matched_sequence']} - {str(r['protein_name'])[:36]} ({pdb}; {r['sam']} + {sub_short}) - {place}"
        (fd / f"{folder}.pdb").write_text("".join(mc.subset(tmp, ev["keep"])))
        (fd / "view.cxc").write_text(mc.cxc(folder, ev["chain"], ev["sub_chains"], ev["kind"], ev["sub_lig"], True, rng, CLASS_HEX[cls], cap))
        (fd / "view.pml").write_text(mc.pml(folder, ev["chain"], ev["sub_chains"], ev["kind"], ev["sub_lig"], True, rng, CLASS_HEX[cls], cap))
        tmp.unlink()
        rows.append({"motif": motif, "consensus": spec.loc[motif, "consensus"], "class": cls, "protein": r["protein_name"], "organism": r["organism"], "accession": r["accession"], "pdb": pdb,
                     "method": method_of(cif), "cofactor": r["sam"], "substrate": sub_short, "chain": ev["chain"], "range": rng, "matched_seq": r["matched_sequence"],
                     "d_SAM": round(ev["d_cof"], 1), "d_substrate": round(ev["d_sub"], 1), "placement": place, "folder": f"{cls}/{folder}",
                     "d_substrate_SAM": round(ev["d_ss"], 1), "resolution": r["resolution"], "in_class": bool(r["in_class"]), "reviewed": bool(r["reviewed"]), "tried": tried,
                     "substrate_kind_fallback": kind != EXPECT[cls], "n_passing": len(passing), "majority_placement": maj})
        print(f"  {motif:13} {pdb} {ev['chain']}:{rng:9} {str(r['matched_sequence']):20} SAM {ev['d_cof']:4.1f}  sub {ev['d_sub']:5.1f}  -> {place}   ({'in-class' if r['in_class'] else 'OUT-OF-CLASS'}, {'reviewed' if r['reviewed'] else 'unreviewed'}, {len(passing)} passing of {tried} tried{', substrate kind FALLBACK: ' + kind if kind != EXPECT[cls] else ''})", flush=True)
    pd.DataFrame(log).to_csv(out / "motif_pick_log.csv", index=False)
    pd.DataFrame(alts).to_csv(out / "motif_pick_alternatives.csv", index=False)
    res = pd.DataFrame(rows)
    res.to_csv(out / "motif_placements.csv", index=False)
    print(f"\n{len(res)} of {cand['motif'].nunique()} motifs have a passing structure -> {out}/motif_placements.csv")
    if len(res):
        print(res["placement"].value_counts().to_string())


if __name__ == "__main__":
    main()
