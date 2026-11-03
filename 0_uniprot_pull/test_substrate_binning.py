#!/usr/bin/env python3
"""Golden cases for substrate_binning.py: real UniProt reaction texts (release 2026_03) with the class a reviewer would give.
Run:  python test_substrate_binning.py   (or pytest)"""
from substrate_binning import bin_entry, classify_participant, label_from_names, parse_reactions, split_participants, transfer_kind

CUR = "Evidence={ECO:0000269|PubMed:1};"
SIM = "Evidence={ECO:0000250|UniProtKB:P12345};"


def ca(reaction, xref="Xref=Rhea:RHEA:1, ChEBI:CHEBI:59789;", ev=CUR):
    return f"CATALYTIC ACTIVITY: Reaction={reaction}; {xref} EC=2.1.1.-; {ev}"


COMP = "Xref=Rhea:RHEA:2, Rhea:RHEA-COMP:1, Rhea:RHEA-COMP:2, ChEBI:CHEBI:59789;"
SAH = "S-adenosyl-L-homocysteine + H(+)"
CASES = [  # (name, catalytic activity text, expected labels)
    ("PRMT arginine", ca(f"L-arginyl-[protein] + S-adenosyl-L-methionine = N(omega)-methyl-L-arginyl-[protein] + {SAH}", COMP), {"Protein"}),
    ("histone lysine", ca(f"L-lysyl-[histone H3] + S-adenosyl-L-methionine = N(6)-methyl-L-lysyl-[histone H3] + {SAH}", COMP), {"Protein"}),
    ("calmodulin (bracket first)", ca(f"[calmodulin]-L-lysine + S-adenosyl-L-methionine = [calmodulin]-N(6)-methyl-L-lysine + {SAH}", COMP), {"Protein"}),
    ("PP2A leucine", ca(f"[phosphatase 2A protein]-C-terminal L-leucine + S-adenosyl-L-methionine = [phosphatase 2A protein]-C-terminal L-leucine methyl ester + S-adenosyl-L-homocysteine", COMP), {"Protein"}),
    ("citrate synthase lysine", ca(f"L-lysyl-[citrate synthase] + S-adenosyl-L-methionine = N(6)-methyl-L-lysyl-[citrate synthase] + {SAH}", COMP), {"Protein"}),
    ("cypemycin precursor peptide", ca(f"N-terminal L-alanyl-[cypemycin] + 2 S-adenosyl-L-methionine = N-terminal N,N-dimethyl-L-alanyl-[cypemycin] + 2 S-adenosyl-L-homocysteine + 3 H(+)", COMP), {"Protein"}),
    ("ACP-tethered", ca(f"malonyl-[ACP] + S-adenosyl-L-methionine = malonyl-[ACP] methyl ester + {SAH}", COMP), {"Protein"}),
    ("tRNA", ca(f"guanosine(37) in tRNA + S-adenosyl-L-methionine = N(1)-methylguanosine(37) in tRNA + {SAH}", COMP), {"RNA"}),
    ("16S rRNA", ca(f"adenosine(1518)/adenosine(1519) in 16S rRNA + 4 S-adenosyl-L-methionine = N(6)-dimethyladenosine in 16S rRNA + 4 {SAH}", COMP), {"RNA"}),
    ("DNA", ca(f"a 2'-deoxycytidine in DNA + S-adenosyl-L-methionine = a 5-methyl-2'-deoxycytidine in DNA + {SAH}", COMP), {"DNA"}),
    ("catechol (COMT)", ca(f"a catechol + S-adenosyl-L-methionine = a guaiacol + {SAH}"), {"Metabolite"}),
    ("no Rhea id (plant OMT)", "CATALYTIC ACTIVITY: Reaction=S-adenosyl-L-methionine + a 3'-hydroxyflavonoid = S-adenosyl-L-homocysteine + a 3'-methoxyflavonoid.; EC=2.1.1.267; " + CUR, {"Metabolite"}),
    ("arsenic: thioredoxin is a redox partner, not the acceptor",
     ca("arsenic triglutathione + [thioredoxin]-dithiol + S-adenosyl-L-methionine + 2 H2O = methylarsonous acid + [thioredoxin]-disulfide + 3 glutathione + S-adenosyl-L-homocysteine + H(+)", COMP), {"Metabolite"}),
    ("sugar polymer is not a protein",
     ca(f"3-O-phospho-alpha-D-Man-(1->2)-[alpha-D-Man-(1->3)-alpha-D-Man-(1->2)](n)-alpha-D-GlcNAc-undecaprenyl diphosphate + S-adenosyl-L-methionine = 3-O-methyl-alpha-D-Man-(1->2)-[alpha-D-Man-(1->3)-alpha-D-Man-(1->2)](n)-alpha-D-GlcNAc-undecaprenyl diphosphate + {SAH}", COMP), {"Metabolite"}),
    ("carboxy-SAM synthase (CmoA): SAM + prephenate", ca("prephenate + S-adenosyl-L-methionine = carboxy-S-adenosyl-L-methionine + 3-phenylpyruvate + H2O"), {"Metabolite"}),
    ("nicotianamine synthase: SAM is the only reactant, so the substrate is a small molecule (the archived pipeline dropped these)", ca("3 S-adenosyl-L-methionine = nicotianamine + 3 S-methyl-5'-thioadenosine + 3 H(+)"), {"Metabolite"}),
    ("two reactions (CARNMT1): metabolite + protein",
     ca(f"carnosine + S-adenosyl-L-methionine = anserine + {SAH}") + " " + ca(f"L-histidyl-[protein] + S-adenosyl-L-methionine = N(pros)-methyl-L-histidyl-[protein] + {SAH}", COMP), {"Metabolite", "Protein"}),
    ("PKS-tethered polyketide", ca(f"(2E,4E)-2,4-dimethyl-6-oxo-dodeca-2,4-dienoyl-[lovastatin nonaketide synthase] + S-adenosyl-L-methionine = (2E,4E)-2,4,6-trimethyl-dodeca-2,4-dienoyl-[lovastatin nonaketide synthase] + {SAH}", COMP), {"Protein"}),
    ("tRNA with a metabolite co-substrate: macromolecule wins within a reaction", ca("L-lysine + tRNA(Lys) + S-adenosyl-L-methionine = L-lysyl-tRNA(Lys) + S-adenosyl-L-homocysteine", COMP), {"RNA"}),
    ("SAM written with a stereo prefix (seen once in UniProt 2026_03): still a SAM reaction, and SAM is not the acceptor",
     ca(f"3-O-acetyl-4'-O-demethylpapaveroxine + (S)-S-adenosyl-L-methionine = 3-O-acetylpapaveroxine + {SAH}"), {"Metabolite"}),
    ("not a SAM reaction (hydrolysis)", "CATALYTIC ACTIVITY: Reaction=Hydrolysis of the 21-Trp-|-Val-22 bond in big endothelin to form endothelin 1.; EC=3.4.24.71; " + CUR, set()),
    ("CmoB: carboxy-SAM is the donor, not SAM", ca("5-hydroxyuridine(34) in tRNA(Ser) + carboxy-S-adenosyl-L-methionine = 5-carboxymethoxyuridine(34) in tRNA(Ser) + S-adenosyl-L-homocysteine + H(+)", "Xref=Rhea:RHEA:3, ChEBI:CHEBI:134278;"), set()),
]


def test_cases():
    bad = []
    for name, text, want in CASES:
        got = bin_entry(text)["labels"]
        if got != want:
            bad.append((name, want, got))
    assert not bad, bad


def test_participants():
    assert classify_participant("an external acceptor")[0] == "Metabolite"  # 'rna' inside a word is not RNA
    assert classify_participant("a nucleoside in ssDNA")[0] == "DNA"
    assert classify_participant("(2S)-2-amino-4-{[(1S)-1-carboxy-2-(1H-imidazol-4-yl)ethyl]amino}butanoate")[0] == "Metabolite"
    assert classify_participant("5'-[hydroxy(phospho)]adenosine", has_comp=False)[0] == "Metabolite"
    assert classify_participant("L-seryl-[protein]", has_comp=False)[0] == "Protein"
    assert split_participants("2 S-adenosyl-L-methionine + n S-adenosyl-L-methionine + a [protein]-L-lysine") == ["S-adenosyl-L-methionine"] * 2 + ["a [protein]-L-lysine"]
    assert split_participants("a (2S)-x(1+2) + H(+)") == ["a (2S)-x(1+2)", "H(+)"]


def test_evidence_and_transfer():
    sim = bin_entry(ca("a catechol + S-adenosyl-L-methionine = a guaiacol + S-adenosyl-L-homocysteine + H(+)", ev=SIM))
    assert sim["source"] == "reaction_uncurated" and sim["labels"] == {"Metabolite"}
    mixed = bin_entry(ca("a catechol + S-adenosyl-L-methionine = a guaiacol + S-adenosyl-L-homocysteine + H(+)") + " " +
                      ca("L-lysyl-[protein] + S-adenosyl-L-methionine = N(6)-methyl-L-lysyl-[protein] + S-adenosyl-L-homocysteine + H(+)", COMP, SIM))
    assert mixed["source"] == "reaction_curated" and mixed["labels"] == {"Metabolite"} and mixed["n_ignored"] == 1  # the gate is the reaction's: the by-similarity protein reaction does not count
    only_other = "CATALYTIC ACTIVITY: Reaction=a thing = another thing; EC=1.1.1.1; " + CUR + " " + ca("a catechol + S-adenosyl-L-methionine = a guaiacol + S-adenosyl-L-homocysteine + H(+)", ev=SIM)
    assert bin_entry(only_other)["source"] == "reaction_uncurated"  # a curated NON-SAM reaction does not vouch for the uncurated SAM one
    assert transfer_kind(parse_reactions(ca("prephenate + S-adenosyl-L-methionine = carboxy-S-adenosyl-L-methionine + 3-phenylpyruvate + H2O"))[0]) == "carboxy"
    assert transfer_kind(parse_reactions(ca("3 S-adenosyl-L-methionine = nicotianamine + 3 S-methyl-5'-thioadenosine + 3 H(+)"))[0]) == "aminopropyl"
    assert bin_entry(None)["source"] == "none" and bin_entry("")["labels"] == set()


def test_all_mode_redox_and_polymers():
    fd = ca("2 reduced [2Fe-2S]-[ferredoxin] + 2 H(+) + a thing = 2 oxidized [2Fe-2S]-[ferredoxin] + a product", "Xref=Rhea:RHEA:5, Rhea:RHEA-COMP:1, Rhea:RHEA-COMP:2;")
    assert bin_entry(fd, "all")["labels"] == {"Metabolite"}  # ferredoxin is a redox partner
    ligase = ca("(deoxyribonucleotide)n-3'-hydroxyl + 5'-phospho-(deoxyribonucleotide)m + NAD(+) = (deoxyribonucleotide)n+m + AMP + beta-nicotinamide D-nucleotide", "Xref=Rhea:RHEA:6;")
    assert bin_entry(ligase, "all")["labels"] == {"DNA"}
    aars = ca("ATP + L-glutamate + tRNA(Glu) = AMP + diphosphate + L-glutamyl-tRNA(Glu)", "Xref=Rhea:RHEA:7;")
    assert bin_entry(aars, "all")["labels"] == {"RNA"}


def test_all_mode_gaps_found_against_the_archive():
    e2 = "CATALYTIC ACTIVITY: Reaction=S-ubiquitinyl-[E1 ubiquitin-activating enzyme]-L-cysteine + [E2 ubiquitin-conjugating enzyme]-L-cysteine = [E1 ubiquitin-activating enzyme]-L-cysteine + S-ubiquitinyl-[E2 ubiquitin-conjugating enzyme]-L-cysteine.; EC=2.3.2.23; " + CUR
    assert bin_entry(e2, "all")["labels"] == {"Protein"}  # no Rhea id: the brackets are trusted
    seca = "CATALYTIC ACTIVITY: Reaction=ATP + H2O + cellular proteinSide 1 = ADP + phosphate + cellular proteinSide 2.; EC=7.4.2.8; Evidence={ECO:0000255|HAMAP-Rule:MF_01382};"
    assert bin_entry(seca, "all")["labels"] == {"Protein"}
    ppk = ca("ATP + [phosphate](n) = ADP + [phosphate](n+1)", "Xref=Rhea:RHEA:8, Rhea:RHEA-COMP:1, Rhea:RHEA-COMP:2;")
    assert bin_entry(ppk, "all")["labels"] == {"Metabolite"}  # polyphosphate is a polymer but not a protein
    iscs = ca("L-cysteine + [sulfur carrier] = L-alanine + [sulfur carrier]-S-sulfanyl", "Xref=Rhea:RHEA:9, Rhea:RHEA-COMP:1, Rhea:RHEA-COMP:2;")
    assert bin_entry(iscs, "all")["labels"] == {"Protein"}  # a protein sulfur carrier is the acceptor, unlike a redox partner


def test_structural_typing_with_rhea_families():
    """Rhea's family decides, so a polymer with a name the word rules cannot know is still typed; the word rules only catch names Rhea lacks."""
    rhea = {"S-adenosyl-L-methionine": "SmallMolecule", "S-adenosyl-L-homocysteine": "SmallMolecule", "H(+)": "SmallMolecule", "Na(+)": "SmallMolecule",
            "L-lysyl-(zeta-polymer 7)": "Polypeptide", "N(6)-methyl-L-lysyl-(zeta-polymer 7)": "Polypeptide", "guanosine(37) in tRNA(Phe)": "Polynucleotide",
            "RNA(n)": "Polymer", "[phosphate](n)": "Polymer", "L-seryl-[protein]": "Polypeptide"}
    new_poly = ca(f"L-lysyl-(zeta-polymer 7) + S-adenosyl-L-methionine = N(6)-methyl-L-lysyl-(zeta-polymer 7) + {SAH}", COMP)
    assert bin_entry(new_poly)["labels"] == {"Metabolite"}  # the word rules alone cannot know it ...
    r = bin_entry(new_poly, "sam", rhea)
    assert r["labels"] == {"Protein"} and r["typed"][0][2].startswith("rhea:Polypeptide")  # ... Rhea's family does
    assert classify_participant("Na(+)(out)", True, rhea) == ("Metabolite", "rhea:SmallMolecule")  # the transport suffix is removed
    assert classify_participant("guanosine(37) in tRNA(Phe)", True, rhea)[0] == "RNA"
    assert classify_participant("RNA(n)", True, rhea)[1].startswith("rhea:Polymer>") and classify_participant("RNA(n)", True, rhea)[0] == "RNA"
    assert classify_participant("[phosphate](n)", True, rhea)[0] == "Metabolite"  # Rhea Polymer family: the word rules decide inside it
    # Rhea files some nucleic-acid residues and a small molecule under its polypeptide class: a nucleic-acid word OUTSIDE a bracket wins, inside one it only names a protein
    rhea.update({"5-methylaminomethyl-2-(Se-phospho)selenouridine(34) in tRNA": "Polypeptide", "5-aminomethyl-dUMP in DNA": "Polypeptide",
                 "L-histidyl-[RNA ligase]": "Polypeptide", "phospho-[DNA-directed RNA polymerase]": "Polypeptide", "a Fe(III)-siderophore": "Polypeptide"})
    assert classify_participant("5-methylaminomethyl-2-(Se-phospho)selenouridine(34) in tRNA", True, rhea)[0] == "RNA"
    assert classify_participant("5-aminomethyl-dUMP in DNA", True, rhea)[0] == "DNA"
    assert classify_participant("L-histidyl-[RNA ligase]", True, rhea)[0] == "Protein"
    assert classify_participant("phospho-[DNA-directed RNA polymerase]", True, rhea)[0] == "Protein"
    assert classify_participant("a Fe(III)-siderophore(in)", True, rhea)[0] == "Metabolite"
    unknown = classify_participant("a thiopurine", False, rhea)  # no Rhea id for this reaction in UniProt: word rules, logged as text:
    assert unknown[0] == "Metabolite" and unknown[1].startswith("text:")


def test_rhea_correcting_a_misfiled_family_changes_nothing():
    """If Rhea moves a mis-filed compound to the right family the label stays; the conflict rule and the override just stop being needed."""
    wrong = {"5-aminomethyl-dUMP in DNA": "Polypeptide", "an N(7)-methylguanosine in RNA": "Polypeptide", "a Fe(III)-siderophore": "Polypeptide"}
    fixed = {"5-aminomethyl-dUMP in DNA": "Polynucleotide", "an N(7)-methylguanosine in RNA": "Polynucleotide", "a Fe(III)-siderophore": "SmallMolecule"}
    want = {"5-aminomethyl-dUMP in DNA": "DNA", "an N(7)-methylguanosine in RNA": "RNA", "a Fe(III)-siderophore": "Metabolite"}
    import substrate_binning as sb
    saved = dict(sb.FAMILY_OVERRIDES)
    try:
        for n, w in want.items():
            assert classify_participant(n, True, wrong)[0] == w            # today: the conflict rule / override catch it
            sb.FAMILY_OVERRIDES.clear()
            assert classify_participant(n, True, fixed)[0] == w            # after Rhea's fix, with no override: same label
            sb.FAMILY_OVERRIDES.update(saved)
            assert classify_participant(n, True, fixed)[0] == w            # and with the override left in place: still the same
    finally:
        sb.FAMILY_OVERRIDES.clear()
        sb.FAMILY_OVERRIDES.update(saved)


def test_carrier_only_flag():
    acp = ca(f"malonyl-[ACP] + S-adenosyl-L-methionine = malonyl-[ACP] methyl ester + {SAH}", COMP)
    assert bin_entry(acp)["carrier_only"] is True
    assert bin_entry(ca(f"L-lysyl-[histone H3] + S-adenosyl-L-methionine = N(6)-methyl-L-lysyl-[histone H3] + {SAH}", COMP))["carrier_only"] is False
    both = ca(f"malonyl-[ACP] + S-adenosyl-L-methionine = malonyl-[ACP] methyl ester + {SAH}", COMP) + " " + ca(f"L-lysyl-[histone H3] + S-adenosyl-L-methionine = N(6)-methyl-L-lysyl-[histone H3] + {SAH}", COMP)
    assert bin_entry(both)["carrier_only"] is False  # one real protein acceptor is enough


def test_sam_name_variants():
    r = bin_entry(ca(f"3-O-acetyl-4'-O-demethylpapaveroxine + (S)-S-adenosyl-L-methionine = 3-O-acetylpapaveroxine + {SAH}"))
    assert [n for n, _, _ in r["typed"]] == ["3-O-acetyl-4'-O-demethylpapaveroxine"] and r["transfer"] == {"methyl"}
    # a different adenosyl compound is NOT SAM: decarboxylated SAM (spermidine synthase) and carboxy-SAM as the donor stay out of the SAM set
    assert bin_entry(ca("S-adenosyl 3-(methylsulfanyl)propylamine + spermidine = thermospermine + S-methyl-5'-thioadenosine + H(+)"))["labels"] == set()
    assert bin_entry(ca("5-hydroxyuridine(34) in tRNA(Ser) + carboxy-S-adenosyl-L-methionine = 5-carboxymethoxyuridine(34) in tRNA(Ser) + S-adenosyl-L-homocysteine + H(+)"))["labels"] == set()


def test_all_mode():
    assert bin_entry("CATALYTIC ACTIVITY: Reaction=Hydrolysis of the 21-Trp-|-Val-22 bond in big endothelin to form endothelin 1.; EC=3.4.24.71; " + CUR, "all")["labels"] == {"Protein"}
    atp = ca("ATP + L-seryl-[protein] = ADP + O-phospho-L-seryl-[protein] + H(+)", "Xref=Rhea:RHEA:4, Rhea:RHEA-COMP:1, Rhea:RHEA-COMP:2, ChEBI:CHEBI:30616;")
    assert bin_entry(atp, "all")["labels"] == {"Protein"}


def test_names():
    assert label_from_names("Ribosomal RNA small subunit methyltransferase G (EC 2.1.1.170)")[0] == {"RNA"}
    assert label_from_names("Protein arginine N-methyltransferase 7")[0] == {"Protein"}
    assert label_from_names("Ribosomal lysine N-methyltransferase 5")[0] == {"Protein"}
    assert label_from_names("Modification methylase HhaI (EC 2.1.1.37)")[0] == {"DNA"}
    assert label_from_names("DNA (cytosine-5)-methyltransferase 3A")[0] == {"DNA"}
    assert label_from_names("Catechol O-methyltransferase")[0] == set()
    assert label_from_names("tRNA/rRNA methyltransferase and histone methyltransferase")[0] == set()  # ambiguous: no label


if __name__ == "__main__":
    for t in (test_cases, test_participants, test_evidence_and_transfer, test_all_mode_redox_and_polymers, test_all_mode_gaps_found_against_the_archive, test_structural_typing_with_rhea_families, test_rhea_correcting_a_misfiled_family_changes_nothing, test_carrier_only_flag, test_sam_name_variants, test_all_mode, test_names):
        t()
        print("ok", t.__name__)
