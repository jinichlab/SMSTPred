#!/usr/bin/env python3
"""For the four motifs with no substrate-bound structure (Protein 2/3/8, Other 2):

  (a) identify the protein whose structure we used for the figure, and
  (b) find the most-characterized protein that carries the motif -- i.e. reviewed
      (SwissProt) FIMO carriers that also have experimental PDB entries.

Uses the real per-protein FIMO scan (fimo_all_{cls}/fimo.tsv) over all class members,
joined to corpus.csv (source=train => reviewed SwissProt), then UniProt REST for
names / organisms / PDB cross-references. Network: login node only.
"""
import io
import json
import urllib.parse
import urllib.request
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
API = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = "accession,reviewed,protein_name,organism_name,length,xref_pdb"

REPS = {"Protein motif 2": "Q6C4Y5", "Protein motif 3": "A7GSD9",
        "Protein motif 8 (NPPY)": "Q9A9T7", "Other motif 2": "(PDB 3C3Y)"}
TARGETS = [("Protein", "2-GGTGLKVWEAAJVLA", "Protein motif 2"),
           ("Protein", "3-AADVVYDP", "Protein motif 3"),
           ("Protein", "8-SNPPYIPT", "Protein motif 8 (NPPY)"),
           ("Other", "2-TGIDLSPEMJEVA", "Other motif 2")]


def uniprot(accs):
    """Fetch UniProt rows for a list of accessions (chunked)."""
    rows = []
    for i in range(0, len(accs), 90):
        chunk = accs[i:i + 90]
        q = " OR ".join(f"accession:{a}" for a in chunk)
        url = f"{API}?{urllib.parse.urlencode({'query': q, 'fields': FIELDS, 'format': 'tsv', 'size': 500})}"
        req = urllib.request.Request(url, headers={"User-Agent": "mtase-motif-check"})
        with urllib.request.urlopen(req, timeout=90) as r:
            txt = r.read().decode()
        if txt.strip():
            rows.append(pd.read_csv(io.StringIO(txt), sep="\t"))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main():
    corpus = pd.read_csv(HERE / "corpus.csv", usecols=["Entry", "source"], low_memory=False)
    src = dict(zip(corpus.Entry, corpus.source))
    fimo = {c: pd.read_csv(HERE / f"fimo_all_{c}" / "fimo.tsv", sep="\t", comment="#")
            .dropna(subset=["sequence_name"]) for c in ("Protein", "Other")}

    # --- (a) the representative proteins we used ---
    reps = [a for a in REPS.values() if a.startswith(("Q", "A", "P", "O"))]
    print("=" * 100)
    print("(a) PROTEINS WHOSE STRUCTURE WE USED")
    ri = uniprot(reps)
    for label, acc in REPS.items():
        if not acc.startswith(("Q", "A", "P", "O")):
            print(f"\n{label}: {acc} -- experimental structure, not an AlphaFold model")
            continue
        r = ri[ri.iloc[:, 0] == acc]
        if len(r):
            r = r.iloc[0]
            pdb = str(r.get("PDB", "") or "").strip()
            print(f"\n{label}: {acc}")
            print(f"   {r['Protein names']}")
            print(f"   organism: {r['Organism']}   length: {r['Length']}")
            print(f"   status: {r['Reviewed']}   PDB: {pdb if pdb else 'none'}")
            print(f"   in our corpus as: {src.get(acc, 'NOT IN CORPUS')}")

    # --- (b) best-characterized carriers ---
    print("\n" + "=" * 100)
    print("(b) MOST-CHARACTERIZED PROTEINS CARRYING EACH MOTIF")
    for cls, mid, label in TARGETS:
        hits = fimo[cls][fimo[cls].motif_id == mid]
        best_p = hits.sort_values("p-value").drop_duplicates("sequence_name")
        carriers = best_p.sequence_name.tolist()
        reviewed = [e for e in carriers if src.get(e) == "train"]
        print(f"\n### {label}  [{mid}]")
        print(f"    carriers {len(carriers):,} | reviewed {len(reviewed)}")
        if not reviewed:
            print("    no reviewed carrier")
            continue
        info = uniprot(reviewed)
        if info.empty:
            print("    (uniprot lookup empty)")
            continue
        info = info.rename(columns={info.columns[0]: "Entry"})
        info["n_pdb"] = info["PDB"].fillna("").apply(lambda s: len([x for x in str(s).split(";") if x.strip()]))
        pmap = dict(zip(best_p.sequence_name, best_p["p-value"]))
        info["fimo_p"] = info.Entry.map(pmap)
        withpdb = info[info.n_pdb > 0].sort_values("n_pdb", ascending=False)
        print(f"    reviewed carriers WITH an experimental PDB: {len(withpdb)}")
        for _, r in withpdb.head(8).iterrows():
            print(f"      {r.Entry}  {r.n_pdb:>3} PDB  p={r.fimo_p:.1e}  {str(r['Protein names'])[:62]}")
            print(f"               {str(r['Organism'])[:60]}  | PDB: {str(r['PDB'])[:70]}")
        if withpdb.empty:
            print("      -- none; top reviewed carriers by FIMO p-value:")
            for _, r in info.sort_values("fimo_p").head(6).iterrows():
                print(f"      {r.Entry}  p={r.fimo_p:.1e}  {str(r['Protein names'])[:60]} | {str(r['Organism'])[:40]}")


if __name__ == "__main__":
    main()
