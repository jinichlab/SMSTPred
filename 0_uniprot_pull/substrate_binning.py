#!/usr/bin/env python3
"""
Stage 0 -- substrate typing ("binning") of UniProt enzymes from their reaction annotation.

Assigns every enzyme the classes of the molecule it acts on: Metabolite / Protein / RNA / DNA (an enzyme can have several).
Pure functions, no I/O and no network (the table reader / writer is bin_substrates.py, the unit tests are test_substrate_binning.py).

Input is the UniProt "Catalytic activity" text of a pull (pull_uniprot.py), one block per reaction:

    CATALYTIC ACTIVITY: Reaction=L-lysyl-[histone] + S-adenosyl-L-methionine = N(6)-methyl-L-lysyl-[histone] + S-adenosyl-L-homocysteine + H(+);
    Xref=Rhea:RHEA:..., Rhea:RHEA-COMP:..., ChEBI:CHEBI:...; EC=2.1.1.-; Evidence={ECO:0000269|PubMed:...};

Rules, applied per reaction and joined per enzyme (a class is set when ANY reaction of the enzyme acts on it):

 1. The reaction is a SAM reaction when S-adenosyl-L-methionine (ChEBI:59789) is a reactant. SAM is a reactant of methyl transfer
    (-> S-adenosyl-L-homocysteine), but also of carboxy-SAM synthesis (CmoA), aminocarboxypropyl transfer (nicotianamine / wybutosine
    synthesis) ...; all of them stay in the set (as in the archived labelling) and `transfer` records which kind it was.
 2. The acceptor(s) are the reactants that are not SAM / water / protons / nucleotide and redox cofactors / redox partner proteins
    (thioredoxin, glutaredoxin, ferredoxin, flavodoxin, glutathione).  [the old pipeline called the arsenite methyltransferases "Protein"
    because [thioredoxin] is in the reaction; the acceptor is arsenic triglutathione]
 3. Each acceptor is typed by Rhea's own classification of the compound (rhea_compounds.tsv.gz from fetch_rhea_compounds.py, looked up by the exact
    name the UniProt text uses; "Na(+)(out)" is looked up as "Na(+)"):
      SmallMolecule   -> Metabolite
      Polypeptide     -> Protein   (a protein-bound residue such as L-histidyl-[protein]; whatever the polymer is called; a new polymer in a later
                                    Rhea release needs no new rule). Redox partners are removed beforehand (see 2).
      Polynucleotide  -> RNA or DNA (the name carries an RNA or DNA word in all 712 compounds of the 2026-11-02 table)
      Polymer         -> the word rules below (Rhea's Polymer family mixes RNA(n), DNA(n), polyglutamylated proteins, sugar chains, polyphosphate, polyprenyls)
    A name that is not in the table (UniProt reactions that have no Rhea id) is typed by the word rules, first match wins, and is logged:
      RNA        an RNA word:   tRNA, rRNA, mRNA, snRNA, "in 16S rRNA", "RNA(n)" ...      (token-bounded; "external" does not match)
      DNA        a DNA word:    "in DNA", ssDNA, dsDNA, (deoxyribonucleotide)n ...
      Protein    a polymer bracket residue-[name] that starts with a letter and is not a sugar chain / polyphosphate / redox partner, or the word
                 protein / peptide / histone
      Metabolite everything else
 4. Per reaction a macromolecule acceptor replaces Metabolite; per enzyme the classes of its reactions are joined. Evidence gate, per REACTION:
    a reaction is "curated" when its own Evidence matches HAMAP / PubMed / ProRule / RuleBase / Ref.N (a curated rule or a literature reference;
    not "by similarity" ECO:0000250, not automatic). The labels of an enzyme come only from its curated reactions (label_source reaction_curated);
    an enzyme none of whose reactions is curated is typed from all of them but kept apart (label_source reaction_uncurated). The archived labelling
    gated per enzyme (any curated block let every block in), which let 88 uncurated SAM reactions and 40 enzymes with only uncurated SAM
    reactions through.
 5. Enzymes without a typable SAM reaction (no reaction text at all, or only non-SAM reactions) can be typed from their names by
    `label_from_names` (label_source name_rule): RNA / DNA / Protein only, never Metabolite. Kept apart: it is a weaker tier whose accuracy
    is measured against the reaction-typed enzymes in bin_substrates.py.

Mode "all" (the whole-Swiss-Prot enzymes set) applies rules 2-4 to every reaction instead of the SAM ones, with all reactants and products
as candidates (a reaction can be "act on a protein" in either direction, e.g. a protease produces peptides).
"""
import re

CLASSES = ("Metabolite", "Protein", "RNA", "DNA")
SAM_NAME = "S-adenosyl-L-methionine"
SAH_NAME = "S-adenosyl-L-homocysteine"
SAM_CHEBI = "CHEBI:59789"
CURATED_RE = re.compile(r"HAMAP|PubMed|ProRule|RuleBase|\{[^}]*Ref\.\d+[^}]*\}", re.I)  # the archived gate (case-insensitive, as then)

# ---------------------------------------------------------------- parsing
_FIELD_AFTER_REACTION = r";\s*(?:Xref|EC|Evidence|PhysiologicalDirection)="
_REACTION = re.compile(r"\s*Reaction=(.*?)(?=" + _FIELD_AFTER_REACTION + r"|;?\s*$)", re.S)
_STOICH = re.compile(r"^(?:\d+|\d*n|m|x|\(n[+-]\d+\))\s+(?=\S)")


def split_participants(side):
    """'2 H(+) + a [protein]' -> ['H(+)', 'a [protein]']; splits on ' + ' outside brackets and parentheses."""
    out, depth, cur, i = [], 0, [], 0
    while i < len(side):
        c = side[i]
        if c in "([{":
            depth += 1
        elif c in ")]}":
            depth = max(0, depth - 1)
        if depth == 0 and side.startswith(" + ", i):
            out.append("".join(cur).strip())
            cur, i = [], i + 3
            continue
        cur.append(c)
        i += 1
    out.append("".join(cur).strip())
    return [_STOICH.sub("", p).strip() for p in out if p.strip()]


def parse_reactions(text):
    """The reaction blocks of a 'Catalytic activity' cell: list of dicts (text, left, right, balanced, xref, ec, evidence, rhea, curated)."""
    if not isinstance(text, str):
        return []
    rows = []
    for block in text.split("CATALYTIC ACTIVITY:")[1:]:
        m = _REACTION.match(block)
        if not m:
            continue
        reaction = m.group(1).strip().rstrip(".").strip()
        balanced = " = " in reaction
        left, right = ([], [])
        if balanced:
            a, b = reaction.split(" = ", 1)
            left, right = split_participants(a), split_participants(b)
        xref = (re.search(r"Xref=([^;]*)", block) or [None, ""])[1]
        ec = (re.search(r"EC=([^;]*)", block) or [None, ""])[1].strip()
        ev = (re.search(r"Evidence=(\{.*?\})\s*(?:;|$)", block, re.S) or [None, ""])[1]
        rows.append({"text": reaction, "left": left, "right": right, "balanced": balanced, "xref": xref, "ec": ec, "evidence": ev,
                     "rhea": re.findall(r"Rhea:RHEA:(\d+)", xref), "has_comp": "RHEA-COMP" in xref, "chebi": re.findall(r"CHEBI:(\d+)", xref),
                     "curated": bool(CURATED_RE.search(block))})  # whole block: a 2nd Evidence= can follow PhysiologicalDirection
    return rows


# SAM as UniProt writes it: plain, or with a stereo prefix ("(S)-S-adenosyl-L-methionine", seen once in release 2026_03)
SAM_RE = re.compile(r"^(?:\((?:S|R|S,S)\)-)?S-adenosyl-L-methionine$")


def is_sam(name):
    return bool(SAM_RE.match(name))


def is_sam_reaction(r):
    """SAM is a reactant (carboxy-SAM as the donor, CmoB, is a different compound and is not a SAM reaction, as in the archived labelling)."""
    return any(is_sam(p) for p in r["left"])


def transfer_kind(r):
    """What SAM does in this reaction: methyl (-> SAH), carboxy (carboxy-SAM synthesis), aminopropyl (-> MTA), other."""
    right = " + ".join(r["right"])
    if SAH_NAME in r["right"]:
        return "methyl"
    if "carboxy-S-adenosyl-L-methionine" in right:
        return "carboxy"
    if "S-methyl-5'-thioadenosine" in right:
        return "aminopropyl"
    return "other"


# ---------------------------------------------------------------- typing one participant
# participants that are never the thing an enzyme acts on
UBIQUITOUS = {SAM_NAME, SAH_NAME, "S-methyl-5'-thioadenosine", "H(+)", "H2O", "ATP", "ADP", "AMP", "GTP", "GDP", "UTP", "CTP", "NAD(+)",
              "NADH", "NADP(+)", "NADPH", "CO2", "O2", "hydrogen peroxide", "phosphate", "diphosphate", "CoA", "AH2", "A", "FAD", "FADH2",
              "FMN", "FMNH2", "reduced FMN", "NH4(+)", "Mg(2+)", "Mn(2+)", "Zn(2+)", "Fe(2+)", "Fe(3+)", "K(+)", "Na(+)", "Ca(2+)", "L-methionine",
              "5'-deoxyadenosine", "glutathione", "glutathione disulfide", "a quinone", "a quinol", "a ubiquinone", "a ubiquinol"}
# electron carriers that appear as residue-[carrier] in redox reactions: they are cofactors of the reaction, not what the enzyme acts on
# (the archived labelling called ~4,900 UP enzymes Protein only because of them: ferredoxin 2,313, thioredoxin 448, plastocyanin 434, flavodoxin 419, ...)
REDOX_PARTNER = re.compile(r"\[[^\]]*(?:redoxin|ferredoxin|flavodoxin|adrenodoxin|plastocyanin|azurin|cytochrome|hemoprotein reductase|electron-transfer flavoprotein|"
                           r"rubredoxin|Dsb[A-Z])[^\]]*\]|^(?:reduced |oxidized )?(?:\[?2Fe-2S\]?|\[?4Fe-4S\]?)|^(?:reduced|oxidized) |\(sulfur carrier\)", re.I)

RNA_RE = re.compile(r"(?<![A-Za-z])(?:[a-z]{1,5}|pre-[a-z]+)?RNAs?(?![A-Za-z])|(?<![A-Za-z])(?:poly)?ribonucleotide\)?[nm]\b|\((?:poly)?ribonucleotide\)[nm]|polyribonucleotide")
DNA_RE = re.compile(r"(?<![A-Za-z])(?:[a-z]{1,3})?DNAs?(?![A-Za-z])|\(deoxyribonucleotide\)[nm]|polydeoxyribonucleotide|\(polydeoxyribonucleotide\)")
# a polymer written as residue-[polymer name]; a bracket in a chemical name begins with a locant or a parenthesis, not a letter
BRACKET_RE = re.compile(r"(?:^|(?<=[\s\-(]))\[([A-Za-z][^\[\]]*)\](?=$|[\s\-()])")
# polymers that are not proteins: sugar chains (a glycosidic linkage "(1->3)"), polyphosphate, glutathione-derived phytochelatin, teichoic acid, ...
SUGAR_POLYMER = re.compile(r"\(\d+->\d+\)|^(?:alpha|beta)-[DL]-|^(?:glycogen|starch|cellulose|chitin|(?:1,\d-)?[a-z]*glucan|mannan|xylan|dextran|"
                           r"peptidoglycan|lipid|glycan|heparan|heparosan|alginate|chondroitin|hyaluron|pectin|lignin|poly\(|polysaccharide|oligosaccharide|"
                           r"phosphate$|polyphosphate|Glu\(-Cys\)|[a-z]*glycerylphospho|teichoic)", re.I)
CARRIER = re.compile(r"\bACP\b|carrier|\bPCP\b|ketide synthase|acid synthase|peptide synthetase|glycine-cleavage complex H protein", re.I)  # carrier-tethered substrates (ACP, PCP, sulfur / lipoyl carriers, PKS / NRPS / FAS)
WORD_PROTEIN = re.compile(r"(?<![A-Za-z])(?i:protein|peptide|histone|polypeptide)s?(?![a-z])")  # 'cellular proteinSide 1': the side is a glued superscript


def macromolecule_brackets(name, has_comp=True):
    """Bracket names of polymers in a participant. Rhea marks a residue-[polymer] reactant with a RHEA-COMP id: for a reaction that Rhea lists
    without one (has_comp False) only a protein-word bracket counts, so a letter-initial bracket inside a chemical name is not taken for a
    polymer. (A reaction with no Rhea id at all, such as the ubiquitin E2 reactions, is trusted: type_reaction passes has_comp=True.)"""
    return [t for t in BRACKET_RE.findall(name) if (has_comp or WORD_PROTEIN.search(t)) and not SUGAR_POLYMER.search(t)
            and not REDOX_PARTNER.search(f"[{t}]")]


_COMPARTMENT = re.compile(r"\((?:in|out)\)$")


def load_rhea_table(path):
    """{compound name: family} from fetch_rhea_compounds.py (Polypeptide | Polynucleotide | Polymer | SmallMolecule). A name that two families share is dropped."""
    import gzip
    table = {}
    with gzip.open(path, "rt") as f:
        next(f)
        for line in f:
            _acc, name, fam = line.rstrip("\n").split("\t")
            table[name] = fam if table.get(name, fam) == fam else "Ambiguous"
    return table


# Rhea files a few small molecules under its generic polypeptide class (the class is its default for a compound with an R group)
FAMILY_OVERRIDES = {"a Fe(III)-siderophore": "SmallMolecule", "a Fe(II)-siderophore": "SmallMolecule"}


def rhea_family(name, rhea):
    if not rhea:
        return None
    bare = _COMPARTMENT.sub("", name)  # transport reactions write Na(+)(out), H(+)(in)
    fam = FAMILY_OVERRIDES.get(bare) or rhea.get(name) or rhea.get(bare)
    return None if fam == "Ambiguous" else fam


def _outside_brackets(name):
    return re.sub(r"\[[^\[\]]*\]", "[]", name)


def _text_rules(name, has_comp):
    """The word rules (RNA word, DNA word, polymer bracket, protein word): used for names Rhea does not know and inside Rhea's mixed Polymer family."""
    if RNA_RE.search(name):
        return "RNA", f"RNA word:{RNA_RE.search(name).group(0)}"
    if DNA_RE.search(name):
        return "DNA", f"DNA word:{DNA_RE.search(name).group(0)}"
    br = macromolecule_brackets(name, has_comp)
    if br:
        tag = "carrier" if any(CARRIER.search(t) for t in br) else "bracket"
        return "Protein", f"{tag}:[{br[0]}]"
    if WORD_PROTEIN.search(name) and not SUGAR_POLYMER.search(name):
        return "Protein", f"word:{WORD_PROTEIN.search(name).group(0)}"
    return "Metabolite", "small molecule / other polymer" if SUGAR_POLYMER.search(name) else "small molecule"


def polymer_token(name):
    """The polymer name inside the bracket of a Polypeptide compound, digits folded (histone H3 -> histone H#), for the review list."""
    m = BRACKET_RE.search(name) or re.search(r"\[([^\[\]]+)\]", name)
    return re.sub(r"\d+", "#", m.group(1)) if m else name[:40]


def classify_participant(name, has_comp=True, rhea=None):
    """-> (class, rule) for one reaction participant name. rule says why (kept in the output for audit): 'rhea:...' structural, 'text:...' word rules."""
    fam = rhea_family(name, rhea)
    if fam == "SmallMolecule":
        return "Metabolite", "rhea:SmallMolecule"
    if fam == "Polynucleotide":
        if RNA_RE.search(name):
            return "RNA", "rhea:Polynucleotide (RNA word)"
        if DNA_RE.search(name):
            return "DNA", "rhea:Polynucleotide (DNA word)"
        return ("DNA" if "deoxy" in name.lower() else "RNA"), "rhea:Polynucleotide (no RNA/DNA word: check)"
    if fam == "Polypeptide":
        out = _outside_brackets(name)  # a nucleic-acid word outside a bracket ("in tRNA", "in DNA") is a nucleic-acid residue Rhea filed as polypeptide;
        if RNA_RE.search(out):         # inside a bracket ("[RNA ligase]", "[DNA-directed RNA polymerase]") it only names a protein
            return "RNA", "rhea:Polypeptide>conflict: RNA word outside a bracket"
        if DNA_RE.search(out):
            return "DNA", "rhea:Polypeptide>conflict: DNA word outside a bracket"
        if REDOX_PARTNER.search(name):
            return "Metabolite", "rhea:Polypeptide (redox partner)"
        return "Protein", f"rhea:Polypeptide:[{polymer_token(name)}]"
    cls, rule = _text_rules(name, has_comp)
    return cls, (f"rhea:Polymer>{rule}" if fam == "Polymer" else f"text:{rule}")


def candidates(r, mode):
    """The participants of a reaction that can be what the enzyme acts on."""
    parts = r["left"] if mode == "sam" else r["left"] + r["right"]
    keep = [p for p in parts if p not in UBIQUITOUS and not is_sam(p) and not REDOX_PARTNER.search(p)]
    seen, out = set(), []
    for p in keep:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def type_reaction(r, mode="sam", rhea=None):
    """-> list of (acceptor name, class, rule) for one parsed reaction."""
    if not r["balanced"]:  # free-text reaction such as 'Hydrolysis of the 21-Trp-|-Val-22 bond in big endothelin'
        if mode == "sam":
            return []
        cls, rule = classify_participant(r["text"], True)  # a free-text reaction is never in the Rhea table: word rules
        if cls == "Metabolite" and re.search(r"\b(?:bond|residue|peptide|protein)s?\b", r["text"], re.I):
            cls, rule = "Protein", "free text: bond / residue"
        return [(r["text"][:80], cls, rule)]
    out = [(p, *classify_participant(p, r["has_comp"] or not r["rhea"], rhea)) for p in candidates(r, mode)]
    if not out:
        out = [("", "Metabolite", "no acceptor left after removing cofactors")]
    return out


# ---------------------------------------------------------------- one enzyme
def bin_entry(catalytic_activity, mode="sam", rhea=None):
    """
    Type one enzyme from its 'Catalytic activity' text.
    mode 'sam': only reactions with S-adenosyl-L-methionine as reactant; mode 'all': every reaction. rhea: load_rhea_table(...) or None (word rules only).
    Per reaction a macromolecule (RNA / DNA / Protein) acceptor replaces Metabolite (an aminoacyl-tRNA synthetase acts on tRNA, not on its amino acid);
    per enzyme the classes of the reactions that are used are joined. Evidence gate per reaction: only curated reactions (their own Evidence matches
    HAMAP / PubMed / ProRule / RuleBase / Ref.N) supply the labels ('reaction_curated'); an enzyme with no curated reaction is typed from all of its
    reactions but marked 'reaction_uncurated' so that it can be kept out of the benchmark.
    -> dict: labels (set of classes), source ('reaction_curated' | 'reaction_uncurated' | 'none'), n_reactions (used), n_ignored (uncurated reactions
             of a curated enzyme that did not count), transfer (set), typed (list of (acceptor, class, rule)), acceptors (the same as
             'name [class: rule]' strings), rhea (list), carrier_only (True when every Protein acceptor is a carrier-tethered substrate).
    """
    rx = [r for r in parse_reactions(catalytic_activity) if (mode == "all" or is_sam_reaction(r))]
    res = {"labels": set(), "source": "none", "n_reactions": 0, "n_ignored": 0, "transfer": set(), "acceptors": [], "typed": [], "rhea": [], "carrier_only": False}
    if not rx:
        return res
    cur = [r for r in rx if r["curated"]]
    use, res["source"] = (cur, "reaction_curated") if cur else (rx, "reaction_uncurated")
    res["n_reactions"], res["n_ignored"] = len(use), len(rx) - len(use)
    for r in use:
        typed = type_reaction(r, mode, rhea)
        classes = {c for _, c, _ in typed}
        if len(classes) > 1:
            classes.discard("Metabolite")
        res["labels"] |= classes
        res["typed"].extend(typed)
        res["acceptors"].extend(f"{name[:70]} [{cls}: {rule}]" for name, cls, rule in typed)
        res["rhea"].extend(r["rhea"])
        if mode == "sam":
            res["transfer"].add(transfer_kind(r))
    res["rhea"] = sorted(set(res["rhea"]), key=int)
    res["acceptors"] = list(dict.fromkeys(res["acceptors"]))
    res["typed"] = list(dict.fromkeys(res["typed"]))
    prot = [n for n, c, _ in res["typed"] if c == "Protein"]
    res["carrier_only"] = bool(prot) and all(CARRIER.search(n) for n in prot)  # Protein only because a carrier protein holds the substrate (ACP, PCP, sulfur carrier ...)
    if not res["labels"]:
        res["source"] = "none"
    return res


# ---------------------------------------------------------------- weaker tier: from the names
NAME_RNA = re.compile(r"(?<![A-Za-z])(?:[a-z]{1,5})?RNA(?![A-Za-z])|(?i:fibrillarin|N\(6\)-adenosine-methyltransferase)|\bRrm\d|\bTrm\d|\bRsm[A-Z]\b|\bRlm[A-Z]\b|\bNSUN\d|\bMETTL(?:3|14|5|16)\b")
NAME_DNA = re.compile(r"(?<![A-Za-z])DNA(?![A-Za-z])|(?i:modification methylase|cytosine-specific methyl|adenine-specific methyl)|\bDNMT\d|\bDnmt\d|\bDam\b|\bDcm\b")
NAME_PROTEIN = re.compile(r"\bprotein[- ](?:arginine|lysine|histidine|glutamine|leucine|L-isoaspartate|N-terminal|[A-Z]-methyl|O-methyl|carboxyl)|histone|\bPRMT\d|"
                          r"ribosomal (?:protein|lysine)|calmodulin|\bSET domain|\bDOT1|\bSETD|\bSMYD|release factor|elongation factor|\bEEF1AKMT|\bN-terminal (?:Xaa|methyl)|"
                          r"\b(?:lysine|arginine|histidine|glutamine) (?:N-)?methyltransferase|\bprotein methyltransferase|\bprotein-lysine", re.I)


def label_from_names(protein_names, ec="", function_cc=""):
    """-> (set of classes, rule) from the protein names; RNA / DNA / Protein only. A name that matches several classes gives no label."""
    s = f"{protein_names or ''} {function_cc or ''}"
    # Protein names are the strong signal; the function text is only used for what the names leave open
    hits = {c for c, rx in (("RNA", NAME_RNA), ("DNA", NAME_DNA), ("Protein", NAME_PROTEIN)) if rx.search(str(protein_names or ""))}
    if not hits:
        hits = {c for c, rx in (("RNA", NAME_RNA), ("DNA", NAME_DNA), ("Protein", NAME_PROTEIN)) if rx.search(s)}
    if len(hits) == 1:
        return hits, "name"
    return set(), "ambiguous" if hits else "no match"
