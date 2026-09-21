#!/usr/bin/env python3
"""Rebuild every specific-motif structure panel from the FIMO carrier evidence.

Replaces the hand-picked (and partly wrong) ranges in make_context_views.py. For
each of the 26 class-specific motifs, motif_top_structures.csv gives the best
in-class carrier whose PDB has BOTH a SAM/SAH cofactor and the substrate bound
(from scan_pdb_ligands.py). This script, per motif:

  1. downloads the structure (mmCIF via gemmi -> works for X-ray and cryo-EM),
  2. LOCATES the motif by searching the FIMO matched_sequence in each protein
     chain (author numbering -- no offset guessing),
  3. identifies the SAM/SAH cofactor and the substrate (nucleic chain / partner
     protein / small-molecule ligand),
  4. measures min heavy-atom distance motif->cofactor and motif->substrate and
     classifies PLACEMENT: near SAM / near substrate / near substrate+SAM / distal,
  5. writes {Class}/{folder}/{folder}.pdb (subset) + view.pml + view.cxc, reusing
     make_context_views styling (grey enzyme, class-colored motif, tan nucleic /
     wheat partner, gold substrate, copper cofactor carbons),
  6. records a row for the README structure table.

Network (RCSB) -> login node only. Output: motif_placements.csv + per-motif folders.
"""
import math
import urllib.request
from pathlib import Path
import pandas as pd
import gemmi
import make_context_views as mc     # reuse subset() / cxc() / pml() / COF_LIGS

S = Path(__file__).resolve().parent
CIF = S / "_carrier_src"
CIF.mkdir(exist_ok=True)
PICKS = pd.read_csv(S / "motif_top_structures.csv")

CLASS_HEX = {"DNA": "6a51a3", "RNA": "cb181d", "Protein": "238b45", "Other": "e6550d"}
COF = set(mc.COF_LIGS) | {"SAM", "SAH", "SFG"}
NEAR = 4.0          # <= this to cofactor/substrate = "near" = direct vdW/H-bond
                    # contact; > this to BOTH cofactor and substrate = "distal"


def fetch(pdb):
    f = CIF / f"{pdb}.cif"
    if not f.exists():
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb}.cif", f)
    return f


def one(res):
    t = gemmi.find_tabulated_residue(res.name)
    return t.one_letter_code.upper() if t and t.is_amino_acid() else ""


def locate(model, matched):
    """Return (chain_id, lo, hi, residues) where the motif matched_sequence sits."""
    for chain in model:
        aa = [r for r in chain if one(r)]
        seq = "".join(one(r) for r in aa)
        i = seq.find(matched)
        if i < 0:                                   # allow 1-mismatch fallback
            for j in range(len(seq) - len(matched) + 1):
                if sum(a != b for a, b in zip(seq[j:j+len(matched)], matched)) <= 1:
                    i = j; break
        if i >= 0:
            res = aa[i:i + len(matched)]
            return chain.name, res[0].seqid.num, res[-1].seqid.num, res
    return None


def atoms(residues):
    return [(a.pos.x, a.pos.y, a.pos.z) for r in residues for a in r]


def mindist(A, B):
    return min((math.dist(p, q) for p in A for q in B), default=float("inf"))


GATE = 28.0          # motif->SAM beyond this = cofactor on a different chain/subunit
                     # (e.g. a ribosome-assembly complex) -> reject, try next carrier
MAX_TRY = 8          # candidate carriers to try per motif before giving up


def evaluate(cls, pdb, matched):
    """Download pdb, locate the motif, identify cofactor + substrate, measure
    distances. Returns a dict of everything needed to render, or None if the motif
    can't be located in this structure."""
    cifpath = fetch(pdb)
    st = gemmi.read_structure(str(cifpath))
    st.setup_entities()
    model = st[0]
    loc = locate(model, matched)
    if not loc:
        return None
    chain, lo, hi, mres = loc
    matoms = atoms(mres)
    method = gemmi.cif.read(str(cifpath)).sole_block().find_value("_exptl.method") or "?"
    method = str(method).strip().strip("'\"").replace("X-RAY DIFFRACTION", "X-ray").replace(
        "ELECTRON MICROSCOPY", "cryo-EM").replace("SOLUTION NMR", "NMR")

    # cofactor: the SAM/SAH copy CLOSEST to the motif (a big complex may hold several)
    cof_at = [(a.pos.x, a.pos.y, a.pos.z) for ch in model for res in ch
              if res.name in COF for a in res]
    keep = {chain}
    sub_chains, sub_lig, kind = [], None, None
    if cls in ("DNA", "RNA"):
        kind = "nucleic"
        nuc = [ch for ch in model if any(
            gemmi.find_tabulated_residue(res.name) and
            gemmi.find_tabulated_residue(res.name).is_nucleic_acid() for res in ch)]
        nuc = sorted(nuc, key=lambda ch: mindist(matoms, atoms(list(ch))))
        sub_chains = [nuc[0].name] if nuc else []
        sub_at = atoms(list(model[sub_chains[0]])) if sub_chains else []
    elif cls == "Protein":
        kind = "partner"
        others = [ch for ch in model if ch.name != chain and any(one(res) for res in ch)]
        others = sorted(others, key=lambda ch: mindist(matoms, atoms([r for r in ch if one(r)])))
        sub_chains = [others[0].name] if others else []
        sub_at = atoms([r for r in model[sub_chains[0]] if one(r)]) if sub_chains else []
    else:
        kind = "ligand"
        sub_lig = str(matched)  # placeholder; overwritten by caller's sub_kind
        sub_at = []
    keep |= set(sub_chains)
    return dict(st=st, model=model, method=method, chain=chain, lo=lo, hi=hi,
                matoms=matoms, cof_at=cof_at, sub_chains=sub_chains, sub_at=sub_at,
                kind=kind, keep=keep)


def main():
    lig = pd.read_csv(Path("/home/nroethler/Code/motif_analysis/motif_structure_ligands.csv"))
    cand = lig[lig.in_class & lig.has_cof & lig.has_sub].copy()
    model_meta = []
    for motif in sorted(PICKS.motif, key=lambda m: (m.split("-")[0], int(m.split("-")[1]))):
        cls = motif.split("-")[0]
        num = motif.split("-")[1]
        rows = cand[cand.motif == motif].sort_values(["reviewed", "p_value"],
                                                     ascending=[False, True])
        chosen = None
        for _, r in rows.head(MAX_TRY).iterrows():
            pdb = str(r.best_pdb).upper()
            ev = evaluate(cls, pdb, str(r.matched_seq))
            if ev is None:
                continue
            # small-molecule substrate: use the real ligand code, re-measure
            if cls == "Other":
                sub_lig = str(r.sub_kind).split(";")[0]
                ev["sub_lig"] = sub_lig
                ev["sub_at"] = [(a.pos.x, a.pos.y, a.pos.z) for ch in ev["model"]
                                for res in ch if res.name == sub_lig for a in res]
                ev["keep"] |= {ch.name for ch in ev["model"] for res in ch
                               if res.name in (COF | {sub_lig})}
            d_cof = mindist(ev["matoms"], ev["cof_at"])
            d_sub = mindist(ev["matoms"], ev["sub_at"])
            if d_cof <= GATE and d_sub <= GATE:     # both ligands on THIS enzyme -> accept
                chosen = (r, pdb, ev, d_cof, d_sub); break
            print(f"    {motif}: {pdb} SAM {d_cof:.0f}A / sub {d_sub:.0f}A (ligand on another chain) -> next")
        if chosen is None:
            print(f"  [SKIP] {motif}: no candidate with cofactor at the enzyme")
            continue
        r, pdb, ev, d_cof, d_sub = chosen
        near_c, near_s = d_cof <= NEAR, d_sub <= NEAR
        place = ("near substrate/SAM" if near_c and near_s else "near SAM" if near_c
                 else "near substrate" if near_s else "distal")

        folder = f"{cls.lower()}_motif{num}_{pdb}"
        d = S / cls / folder
        d.mkdir(parents=True, exist_ok=True)
        ev["st"].remove_alternative_conformations()
        tmp = d / f"_full_{pdb}.pdb"
        ev["st"].write_pdb(str(tmp))
        rng = f"{ev['lo']}-{ev['hi']}"
        hx = CLASS_HEX[cls]
        sub_lig = ev.get("sub_lig")
        sub_short = {"DNA": "DNA", "RNA": "RNA", "Protein": "partner protein"}.get(cls, sub_lig)
        cap = (f"{motif} {r.matched_seq} - {str(r.protein_name)[:36]} ({pdb}; {ev['method']}; "
               f"{r.cof} + {sub_short}) - {place}")
        (d / f"{folder}.pdb").write_text("".join(mc.subset(tmp, ev["keep"])))
        (d / "view.cxc").write_text(mc.cxc(folder, ev["chain"], ev["sub_chains"], ev["kind"],
                                           sub_lig, True, rng, hx, cap))
        (d / "view.pml").write_text(mc.pml(folder, ev["chain"], ev["sub_chains"], ev["kind"],
                                           sub_lig, True, rng, hx, cap))
        tmp.unlink()

        model_meta.append({
            "motif": motif, "consensus": r.consensus, "class": cls,
            "protein": r.protein_name, "organism": r.organism, "accession": r.accession,
            "pdb": pdb, "method": ev["method"], "cofactor": r.cof, "substrate": sub_short,
            "chain": ev["chain"], "range": rng, "matched_seq": r.matched_seq,
            "d_SAM": round(d_cof, 1), "d_substrate": round(d_sub, 1), "placement": place,
            "folder": f"{cls}/{folder}",
        })
        print(f"  {motif:10} {pdb} {ev['chain']}:{rng:9} {str(r.matched_seq):18} "
              f"SAM {d_cof:4.1f}  sub {d_sub:5.1f}  -> {place}")

    out = pd.DataFrame(model_meta)
    out.to_csv(S / "motif_placements.csv", index=False)
    print(f"\n[wrote] motif_placements.csv ({len(out)}/{len(PICKS)} motifs)")
    print(out.placement.value_counts().to_string())


if __name__ == "__main__":
    main()
