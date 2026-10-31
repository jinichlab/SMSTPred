#!/usr/bin/env python3
"""Evidence-based structure picks for every class-specific STREME motif.

For each of the 26 motifs that pass the specificity filter (own-class rate >= 2x
pooled background, from the FIMO cross-scan), list the FIMO carriers that have an
experimental PDB structure -- so the Figure 9/10 structure panels are chosen from
data, not hand-picked. Replaces the eyeballed ranges in make_context_views.py.

Selection (per user): rank carriers by FIMO p-value regardless of class, but flag
in-class vs out-of-class and reviewed (SwissProt) vs not, because the figures should
use IN-CLASS carriers with a substrate/SAM-bound structure -- with the one exception
worth noting where a motif has NO in-class structure but does in another class.

Pipeline:
  fimo_cross/{cls}.tsv  -> carriers of each motif (best p-value hit per protein,
                          with UniProt start/stop + matched_sequence)
  corpus.csv            -> each carrier's own class + reviewed status (source=train)
  UniProt REST          -> protein name, organism, and xref_pdb (login node only)

Output:
  motif_structure_carriers.csv   one row per (motif, structure-backed carrier)
  motif_structure_carriers.txt   readable per-motif top picks + coverage summary

Validation: the 11 structure panels we already trust should re-surface here as top
in-class carriers of their motif; the mis-mapped ones (e.g. rna_motif2/METTL16)
should NOT appear as carriers at all.
"""
import io
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
API = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = "accession,reviewed,protein_name,organism_name,length,xref_pdb"
CLASSES = ["DNA", "RNA", "Protein", "Metabolite"]
SPEC_MIN = 2.0
# how many carriers per motif to look up (bounds network; structures are enriched
# in reviewed + strong-p hits, so we over-sample those)
TOP_INCLASS = 120
TOP_OUTCLASS = 60
PDB_RE = re.compile(r"\b[0-9][A-Za-z0-9]{3}\b")


def uniprot(accs):
    """UniProt rows for a list of accessions (chunked search; login node only)."""
    rows = []
    accs = sorted(set(accs))
    for i in range(0, len(accs), 90):
        chunk = accs[i:i + 90]
        q = " OR ".join(f"accession:{a}" for a in chunk)
        url = f"{API}?{urllib.parse.urlencode({'query': q, 'fields': FIELDS, 'format': 'tsv', 'size': 500})}"
        for attempt in range(4):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "mtase-motif-structs"})
                with urllib.request.urlopen(req, timeout=120) as r:
                    txt = r.read().decode()
                break
            except Exception as e:
                if attempt == 3:
                    print(f"  [warn] chunk {i} failed: {e}"); txt = ""
                else:
                    time.sleep(3 * (attempt + 1))
        if txt.strip():
            rows.append(pd.read_csv(io.StringIO(txt), sep="\t"))
        if (i // 90) % 10 == 0:
            print(f"  uniprot {i}/{len(accs)}", flush=True)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def load_carriers():
    """Best-p-value hit per (motif, protein) across the whole cross-scan."""
    corpus = pd.read_csv(HERE / "corpus.csv", usecols=["Entry", "class", "source"],
                         keep_default_na=False, low_memory=False)
    corpus["class"] = corpus["class"].replace({"None": "Unassigned"})
    prot_class = dict(zip(corpus["Entry"], corpus["class"]))
    reviewed = dict(zip(corpus["Entry"], corpus["source"] == "train"))

    frames = []
    for cls in CLASSES:
        d = pd.read_csv(HERE / "fimo_cross" / f"{cls}.tsv", sep="\t")
        d = d.dropna(subset=["sequence_name"])
        d["motif_class"] = cls
        d["motif"] = cls + "-" + d["motif_id"].astype(str).str.split("-").str[0]
        d["consensus"] = d["motif_id"].astype(str).str.split("-").str[1]
        frames.append(d[["motif", "motif_class", "consensus", "sequence_name",
                         "start", "stop", "p-value", "matched_sequence"]])
    hits = pd.concat(frames, ignore_index=True)
    hits = (hits.sort_values("p-value")
                .drop_duplicates(["motif", "sequence_name"]))     # best hit per protein
    hits["prot_class"] = hits["sequence_name"].map(prot_class)
    hits["in_class"] = hits["prot_class"] == hits["motif_class"]
    hits["reviewed"] = hits["sequence_name"].map(reviewed).fillna(False)
    return hits


def main():
    spec = pd.read_csv(HERE / "motif_specificity_matrix.csv")
    specific = set(spec.loc[spec.ratio_vs_bg >= SPEC_MIN, "motif"])
    print(f"[carriers] {len(specific)} class-specific motifs (>= {SPEC_MIN}x background)")

    hits = load_carriers()
    hits = hits[hits.motif.isin(specific)]

    # pick which carriers to look up: all reviewed in-class + top-p in-class + top-p out
    want = []
    for motif, g in hits.groupby("motif"):
        inc = g[g.in_class]
        out = g[~g.in_class]
        want += list(inc[inc.reviewed].sequence_name)
        want += list(inc.nsmallest(TOP_INCLASS, "p-value").sequence_name)
        want += list(out.nsmallest(TOP_OUTCLASS, "p-value").sequence_name)
    want = sorted(set(want))
    print(f"[carriers] looking up {len(want):,} candidate accessions on UniProt ...")
    meta = uniprot(want)
    if meta.empty:
        print("  [fatal] no UniProt rows returned"); return
    meta.columns = ["accession", "reviewed_up", "protein_name", "organism", "length", "pdb_raw"]
    meta["pdb_raw"] = meta["pdb_raw"].fillna("")
    meta["pdbs"] = meta["pdb_raw"].map(lambda s: sorted(set(PDB_RE.findall(str(s)))))
    meta["n_pdb"] = meta["pdbs"].map(len)
    has_pdb = meta[meta.n_pdb > 0].set_index("accession")

    # build the structure-backed carrier table
    rows = []
    for motif, g in hits.groupby("motif"):
        cons = g.consensus.iloc[0]
        for _, h in g.iterrows():
            acc = h.sequence_name
            if acc not in has_pdb.index:
                continue
            m = has_pdb.loc[acc]
            rows.append({
                "motif": motif, "consensus": cons, "motif_class": h.motif_class,
                "accession": acc, "prot_class": h.prot_class, "in_class": h.in_class,
                "reviewed": bool(h.reviewed), "protein_name": m.protein_name,
                "organism": m.organism, "length": m.length,
                "n_pdb": int(m.n_pdb), "pdbs": ";".join(m.pdbs),
                "fimo_start": int(h.start), "fimo_stop": int(h.stop),
                "matched_seq": h.matched_sequence, "p_value": h["p-value"],
            })
    tab = pd.DataFrame(rows)
    # rank so the best figure candidate per motif floats to the top
    tab = tab.sort_values(
        ["motif", "in_class", "reviewed", "p_value"],
        ascending=[True, False, False, True]).reset_index(drop=True)
    tab.to_csv(HERE / "motif_structure_carriers.csv", index=False)

    # ---- report ----
    L = [f"Structure-backed carriers for {len(specific)} class-specific motifs "
         f"({len(tab)} carrier-structures)\n"]
    for motif in sorted(specific, key=lambda m: (m.split("-")[0], int(m.split("-")[1]))):
        t = tab[tab.motif == motif]
        inc = t[t.in_class]
        L.append(f"{motif} ({tab[tab.motif==motif].consensus.iloc[0] if len(t) else ''})"
                 f"  in-class w/ structure: {len(inc)}   any-class: {len(t)}")
        if len(inc) == 0 and len(t):
            best = t.iloc[0]
            L.append(f"    !! NO in-class structure. best other-class: {best.accession} "
                     f"[{best.prot_class}] {str(best.protein_name)[:45]} PDB {best.pdbs.split(';')[0]}")
        for _, r in inc.head(3).iterrows():
            star = "*" if r.reviewed else " "
            L.append(f"   {star}{r.accession}  {str(r.protein_name)[:46]:46} "
                     f"{r.fimo_start}-{r.fimo_stop} {r.matched_seq:18} PDB {r.pdbs.split(';')[0]} "
                     f"(+{r.n_pdb-1})" if r.n_pdb > 1 else
                     f"   {star}{r.accession}  {str(r.protein_name)[:46]:46} "
                     f"{r.fimo_start}-{r.fimo_stop} {r.matched_seq:18} PDB {r.pdbs.split(';')[0]}")
        L.append("")
    n_missing = sum(1 for m in specific if not len(tab[(tab.motif == m) & tab.in_class]))
    L.append(f"SUMMARY: {len(specific)-n_missing}/{len(specific)} motifs have an in-class "
             f"structure-backed carrier; {n_missing} have none in-class.")
    txt = "\n".join(L)
    (HERE / "motif_structure_carriers.txt").write_text(txt + "\n")
    print(txt[-2500:])
    print("\n[wrote] motif_structure_carriers.csv + .txt")


if __name__ == "__main__":
    main()
