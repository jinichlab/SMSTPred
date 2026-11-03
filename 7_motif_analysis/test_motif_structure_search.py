#!/usr/bin/env python3
"""classify(): the tiers on made-up RCSB records shaped like the real ones (1MHT: M.HhaI + DNA + SAH; an apo entry; a peptide complex; a SAM-only entry)."""
from motif_structure_search import classify


def entry(polymers, ligands):
    return {"exptl": [{"method": "X-RAY DIFFRACTION"}], "rcsb_entry_info": {"resolution_combined": [2.0]}, "struct": {"title": "t"},
            "polymer_entities": [{"entity_poly": {"rcsb_entity_polymer_type": t, "rcsb_sample_sequence_length": n},
                                  "rcsb_polymer_entity_container_identifiers": {"auth_asym_ids": ch, "reference_sequence_identifiers": ([{"database_name": "UniProt", "database_accession": u}] if u else None)}}
                                 for t, n, ch, u in polymers],
            "nonpolymer_entities": [{"nonpolymer_comp": {"chem_comp": {"id": c, "formula_weight": 250.0}}, "rcsb_nonpolymer_entity_container_identifiers": {"auth_asym_ids": ch}} for c, ch in ligands]}


P = "P05102"
hha = entry([("DNA", 12, ["B"], None), ("DNA", 13, ["C"], None), ("Protein", 327, ["A"], P)], [("SAH", ["A"])])
r = classify(hha, P)
assert r["tier"] == "sam_substrate_on_chain" and r["substrate_type"] == "nucleic_acid" and r["sam"] == "SAH", r
assert classify(entry([("Protein", 300, ["A"], P)], [("GOL", ["A"]), ("MG", ["A"])]), P)["tier"] == "apo"                 # buffer / ion only
assert classify(entry([("Protein", 300, ["A"], P)], [("SAM", ["A"])]), P)["tier"] == "sam_only"
assert classify(entry([("Protein", 300, ["A"], P), ("Protein", 12, ["P"], None)], [("SFG", ["A"])]), P)["tier"] == "sam_substrate_on_chain"   # peptide substrate
assert classify(entry([("Protein", 300, ["A"], P), ("Protein", 300, ["B"], None)], [("SAH", ["A"])]), P)["tier"] == "sam_only"            # a long partner chain is not a peptide substrate
small = classify(entry([("Protein", 300, ["A"], P)], [("SAH", ["A"]), ("XYZ", ["A"])]), P)
assert small["tier"] == "sam_substrate_on_chain" and small["substrate"] == "XYZ" and small["substrate_type"] == "small_molecule"
# SAM on another protein of the complex (chain B holds SAH, the carrier is chain A): in the entry, not on the carrier's chain
other = classify(entry([("Protein", 300, ["A"], P), ("Protein", 300, ["B"], "Q00000"), ("DNA", 12, ["C"], None)], [("SAH", ["B"])]), P)
assert other["tier"] == "sam_substrate", other
assert classify(entry([("Protein", 300, ["A"], P)], [("XYZ", ["A"])]), P)["tier"] == "substrate_only"
for junk in ("FES", "UNX", "ATP"):      # iron-sulfur cluster, unknown atom, nucleotide: none of them is a substrate
    assert classify(entry([("Protein", 300, ["A"], P)], [("SAH", ["A"]), (junk, ["A"])]), P)["tier"] == "sam_only", junk
tiny = entry([("Protein", 300, ["A"], P)], [("SAH", ["A"]), ("QQQ", ["A"])])
tiny["nonpolymer_entities"][1]["nonpolymer_comp"]["chem_comp"]["formula_weight"] = 60.0      # too small to be a substrate
assert classify(tiny, P)["tier"] == "sam_only"
print("motif_structure_search classify OK")
