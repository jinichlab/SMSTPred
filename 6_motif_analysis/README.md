# Sequence-Motif, Attribution & Structural Analysis of SAM-MTase Substrate Specificity

This directory contains the **motif / interpretability** stage of the project: a set of
analyses that ask **what sequence and structural features the ProstT5→MLP classifier
uses to distinguish methyltransferase (MTase) substrate type** (DNA / RNA / Protein /
Other-small-molecule). It combines three complementary lines of evidence:

1. **Discriminative linear-motif discovery** (STREME) — what short sequence motifs
   separate each substrate class;
2. **Per-residue attribution** (Integrated Gradients on the classifier) — which residues
   actually drive each call, especially for proteins that have *no* discriminative motif;
3. **Structural mapping** — placing both the motifs and the IG-important residues onto
   experimental (and where necessary AlphaFold / homologous) structures.

Together they converge on a single conclusion (Section 9): **the classifier reads
substrate identity from family-specific buried fold architecture, not from the conserved
catalytic/cofactor machinery.**

---

## 0. Input

The labelled **corpus** of single-substrate SAM-MTases, assembled by `1_build_corpus.py`
from two sources:

| Source | Description | n |
|---|---|---|
| training | SwissProt-reviewed SAM_swiss enzymes (ground-truth labels) | 7,405 |
| deployment | reference-proteome SAMs labelled by the deployment classifier (`call_*`) | ~111k |

`single_label()` keeps proteins with exactly one substrate (Mixed/None dropped) →
**`corpus_all.fasta`** (106,608 sequences) and per-class FASTAs. Every member carries the
SAM-MTase Rossmann superfamily domain (InterPro **IPR029063**) — i.e. all are genuine
methyltransferases. Outputs: `corpus.csv` (Entry, Sequence, InterPro, class, source,
in_domain).

| Class | single-label corpus |
|---|---|
| Other | 41,941 |
| RNA | 39,079 |
| Protein | 18,766 |
| DNA | 6,822 |

---

## 1. Redundancy reduction → family representatives

Run by `run_dedup.sh` (→ `2_make_reps.py`): `mmseqs easy-cluster` at **50 % identity**
collapses paralogues so motif discovery sees one representative per ~family rather than
thousands of near-duplicates.

| Class | family reps |
|---|---|
| Other | 10,470 |
| RNA | 5,986 |
| Protein | 3,959 |
| DNA | 2,561 |
| **Total** | **22,976** |

`reps_<class>.clean.fasta` = the reps with non-standard residues mapped to X (STREME/FIMO
input).

---

## 2. Domain enrichment (orthogonal validation)

`3_domain_enrichment.py` — per-class **InterPro** enrichment (Fisher exact) →
`domain_enrichment.csv`. Confirms each predicted class recovers the expected enzyme
families: **DNA** → C5-cytosine/DNMT + N6-adenine/N4 MTases; **Protein** → PRMT/SET +
HemK; **RNA** → RsmB/NOP2/FtsJ; **Other** → O-MTase / SABATH. This is independent of the
classifier and validates the labels.

---

## 3. Discriminative motif discovery (STREME)

`run_streme.sh` runs **STREME** (MEME suite) on the **full family-rep set** per class
(positives = that class's reps; control = the other three classes' reps), width **8–18**,
**8 motifs/class**. Discriminative = enriched in the class *vs* the others. Outputs:
`streme_<class>/streme.{txt,html}`.

### The 32 motifs — literature status and structural placement

For every motif we record: the consensus and E-value; whether it corresponds to a
**named element in the methyltransferase literature** (the Schubert *et al.* 2003 Motif
I–X scheme, the catalytic-motif and target-recognition-domain literature) or is
**novel/unnamed**; the structure it was mapped onto (Section 8); and **where it physically
sits** — the measured minimum distance from the motif to the **SAM/SAH cofactor** and to
the **substrate** (<4.5 Å = direct contact). "Placement" summarises this as near-SAM,
near-substrate, both (active-site/junction), or neither (a distal structural element).

**DNA** — mapped on DNA-bound crystals (HhaI 1MHT, TaqI 1G38); all motifs have both
cofactor and DNA present, so placement is fully resolved.

| # | Consensus | Literature | Placement (→SAM / →DNA, Å) |
|---|---|---|---|
| 1 | IGG**PCQ**GFSRAG | **KNOWN** — Motif IV catalytic PCQ (C5-cytosine nucleophile, HhaI Cys81) | **active site** — base 1.8, SAM 3.6 |
| 2 | PFA**GSGG**T | **KNOWN** — Motif I AdoMet glycine loop | **near SAM** — 2.9 (DNA 12.7) |
| 3 | YGV**PQRRER**VF | **KNOWN** — target-recognition domain (DNA-contacting basic loop) | **near DNA** — 2.8 (SAM 8.0) |
| 4 | LTPREAARLQGFPDDYTF | novel/unnamed (within the TRD/C-term) | **near SAM** — 3.8 (DNA 9.6, *not* TRD) |
| 5 | KQIGNAVPPPVARA | novel/unnamed | **junction** — SAM 3.2, base 3.4 |
| 6 | FL**ENV**PGL | **KNOWN** — Motif VI catalytic Glu (HhaI Glu119) | **near DNA** — base 2.8 (SAM 5.0) |
| 7 | YTB**PPY**NT | **KNOWN** — Motif IV catalytic [D/N]PPY (N6-adenine, TaqI Asn105) | **active site** — DNA 2.9, SAM 2.9 |
| 8 | GZFFTPREI | novel/unnamed (N-terminal, near cofactor) | **junction** — DNA 3.3, SAM 3.8 |

**Protein** — 2NXN (PrmA + L11 substrate, **no cofactor**), 4GQB (PRMT5 + MEP50 partner +
SAM-analog), AlphaFold for 2/3/5/8 (no ligands).

| # | Consensus | Literature | Placement |
|---|---|---|---|
| 1 | EL**GAGTG**JLGJAAAKL | **KNOWN** — Motif I AdoMet glycine loop (seven-β-strand class) | SAM loop, but **2NXN is apo of cofactor**; 8.2 Å from L11 |
| 2 | **GGTG**LKVWEAAJVLA | **KNOWN** — Motif I / post-I | AlphaFold (no ligand context) |
| 3 | AADVVYDP | tentative — Motif II acidic | AlphaFold |
| 4 | VDIJVSEWMGYFLLNES | novel/unnamed | **active site** — SAM 2.8, partner 2.8 |
| 5 | **THW**KQTVFYL | **KNOWN** — PRMT "THW" loop | AlphaFold |
| 6 | VLIPRPETE | novel/unnamed | **near substrate/partner** — 2.7 (SAM 54, other domain) |
| 7 | AF**GTGT**HPTTRLCLEAL | novel/unnamed (SAM-adjacent) | **near substrate** — L11 3.4 (PrmA apo of SAM) |
| 8 | S**NPPY**IPT | **KNOWN** — glutamine N-MTase catalytic NPPY (HemK/PrmC) | AlphaFold |

**RNA** — 2ZZM (aTrm5 + tRNA + SAM), 9U48 (METTL16 + RNA + SAM), 3M6U (RsmF, **apo**),
AlphaFold for motif 8.

| # | Consensus | Literature | Placement |
|---|---|---|---|
| 1 | DRVLLD**PC**SGLGT | **KNOWN** — m5C catalytic Cys (PC) | RsmF, **apo** (no ligand context) |
| 2 | RJVYS**TC**S | **KNOWN** — m5C catalytic Cys region | distal (METTL16: SAM 11, RNA 19) |
| 3 | RVLDLCA**APGG**KTQV | **KNOWN** — AdoMet-binding region | **near SAM** — 3.0 (RNA 6.4) |
| 4 | LQDASSMLPALALDPQP | novel/unnamed | RsmF, **apo** |
| 5 | TGGFFVAKL | novel/unnamed | RsmF, **apo** |
| 6 | ENEAVVAAALK | novel/unnamed | distal (RNA 11.5; SAM in other domain) |
| 7 | **KKKKKKRK** | functional element — **basic RNA-binding patch** (not a named catalytic motif) | **near RNA** — tRNA 5.2 |
| 8 | **RSR**AAFKL | novel/unnamed (RlmE basic) | AlphaFold |

**Other (small molecule)** — all homology-mapped onto SAH-bound plant O-MTases (1KYZ,
1M6E, 3C3Y); placement is therefore measured on the homolog.

| # | Consensus | Literature | Placement (→SAH / →substrate, Å) |
|---|---|---|---|
| 1 | RVLKPGGRLVL | novel/unnamed | distal — SAH 9.0, sub 14.8 |
| 2 | TGIDLSPEMJEVA | novel/unnamed | **near SAH** — 2.7 (3C3Y has no substrate) |
| 3 | LVDV**GGGTG**QLL | **KNOWN** — Motif I AdoMet glycine loop | **near SAH** — 3.1 (sub 12.2) |
| 4 | GDAEALPFPDGSFDLVVS | novel/unnamed | **near SAH** — 2.8 (sub 9.8) |
| 5 | GGKERTEEEWRAL | novel/unnamed | distal — substrate 6.6, SAH 11.8 |
| 6 | LHDWPDEDCVKIL | novel/unnamed | **active site** — SAH 3.4, substrate 3.1 |
| 7 | HDFFTPQPVKGADVYL | novel/unnamed | **near SAH** — 2.7 (sub 10.7) |
| 8 | GIVQDLPEVIA | novel/unnamed | **near SAH** — 2.5 (sub 11.4) |

**Summary of literature status.** ~**13 of 32** correspond to a named element — the AdoMet
Motif I glycine loop (DNA-2, Prot-1/2, RNA-3, Other-3), the catalytic Motif IV
(DNA-1 PCQ, DNA-7 DPPY, Prot-8 NPPY, RNA-1/2 m5C-Cys), Motif VI catalytic Glu (DNA-6),
the DNA target-recognition domain (DNA-3), the PRMT THW loop (Prot-5), and the basic
RNA-binding patch (RNA-7). The remaining **~19** are not named motifs; structurally most
fall into two camps — **near the SAM/SAH pocket** (a conserved-fold region without a
canonical name, e.g. Other-2/4/7/8) or **distal/structural** (e.g. DNA-4/5, RNA-6,
Other-1/5). Importantly the strongest motif in **every** class is a recognised
catalytic/cofactor signature, validating that STREME found real biology; and putting the
candidates on structure *upgraded* several (e.g. DNA-6 = the catalytic Glu, DNA-8 = a
cofactor-junction element) while showing others are genuine SAM-pocket-adjacent fold
features rather than artefacts.

---

## 4. Motif prevalence (FIMO) — who has a motif, who doesn't

`fimo_<class>/` — each class's 8 motifs scanned back against its reps (FIMO, p<1e-4).
A rep is "no-motif" if it matches **none** of its class's motifs.

| Class | reps | match ≥1 motif | **no motif** | % no motif |
|---|---|---|---|---|
| RNA | 5,986 | 4,367 | **1,619** | **27.0 %** |
| DNA | 2,561 | 2,145 | 416 | 16.2 % |
| Protein | 3,959 | 3,560 | 399 | 10.1 % |
| Other | 10,470 | 9,881 | 589 | 5.6 % |
| **Total** | **22,976** | 19,953 | **3,023** | **13.2 %** |

RNA is the outlier (its motifs are weak/binding-patch rather than a sharp catalytic
signature). These **3,023 no-motif reps** are the focus of the attribution analysis.

---

## 5. Are the no-motif proteins real? (yes)

Characterisation of the 3,023 no-motif reps vs the motif-having majority:

- **100 % carry the SAM-MTase fold** (IPR029063) → they are *genuine methyltransferases*;
  "no motif" ≠ "not real".
- **Enriched for out-of-domain** (27 % in-domain vs 65 % for motif-having) — divergent,
  extrapolated calls; DNA no-motif is 99 % out-of-domain.
- **Enriched for short/fragment sequences** (15 % <150 aa vs 1 %).

So the no-motif set = real MTases whose discriminative signal is *distributed/structural*
rather than a short linear motif (plus a fragment tail).

---

## 6. Per-residue attribution (Integrated Gradients)

`attribution/integrated_gradients.py` + `no_motif_attribution.py`.

**Method.** The classifier eats a **mean-pooled** ProstT5 embedding, so naive saliency is
uniform/useless. IG attributes through each residue's own embedding:
`A_i = (1/L)(e_i − b)·G`, where `G` is the path-averaged gradient of the (logit-space)
class score and `b` is the mean-MTase baseline. Completeness verified (Σ A_i = F(z)−F(b),
ratio ≈ 1.0). **Logit space** is essential — the deployed probability saturates at ~0.999
(vanishing gradient) and the mean baseline is degenerate for the majority class (RNA).

**No-motif aggregate trends** (`no_motif_attr_output/`): per class, the amino-acid
identities over-represented among the top-10 % IG residues, and the peakedness (Gini) of
the attribution:

| Class | top-enriched residues | median Gini |
|---|---|---|
| DNA | **F, T, W** (aromatic) | 0.556 (most distributed) |
| Protein | **L, M, Y** | 0.636 |
| RNA | **N, C, K** (Asn / catalytic Cys / basic) | 0.692 (most focal) |
| Other | **T, W, F** (Thr + aromatic) | 0.613 |

Each class shows a chemically-sensible, *distinct* fingerprint even without a linear
motif. `top_windows_<class>.fasta` holds the ±4 window around each protein's single
most-essential residue (ready to re-STREME for short sub-motifs).

---

## 7. IG on structure — the central result

For one representative **no-motif** protein per class that has an **experimental
structure**, full per-residue IG was mapped onto the structure (B-factor column) and the
burial of the top-10 IG residues quantified (CA neighbours within 10 Å — high = buried
core, low = surface). Deliverables in `structures/IG_<class>/` (`*_IG.pdb` + `view.pml`/
`view.cxc`).

| Class | Protein / structure (cofactor, substrate) | top-10 IG burial | protein avg | substrate-contacting |
|---|---|---|---|---|
| DNA | T4 Dam / **1YF3** (AdoHcy + DNA) | **18.3** | 16.6 | 15.3 (surface) |
| RNA | Q9V2G1 / **5WT1** (SAM + tRNA, cognate) | **20.4** | 17.2 | 14.6 (surface) |
| Protein | NDUFAF7 / **5ZTZ** (apo) | **20.5** | 17.1 | — |
| Other | PsiM / **8PB3** (SAH + norbaeocystin, catalytic 3.5 Å) | **18.8** | 17.7 | distal to both |

**Finding (all 4 classes):** the IG-important residues are **buried hydrophobic core**
residues (coordination 18.3–20.5, above the ~17 average), markedly **more buried than the
substrate-binding residues** (which are surface, ~15, and polar). This holds even in
**cognate/productive complexes** (RNA: SAM–tRNA 4.4 Å; Other: SAH–substrate 3.5 Å), so it
is not an artefact of non-productive crystal forms.

**Interpretation.** All SAM-MTases share the catalytic SAM-binding machinery, so the
active site cannot discriminate substrate type. ProstT5 is *structure-aware* (encodes 3Di
structural tokens); the classifier therefore reads substrate identity from
**family-specific fold architecture encoded in the buried core**, not from the conserved
catalytic or substrate-contacting residues. This explains why the signal is distributed,
why STREME misses many proteins (the signal is structural, not a linear motif), and why
IG repeatedly lands on hydrophobic core clusters.

---

## 8. Motif-on-structure — the structural examples in detail

`structures/<class>/` — one folder per motif, each with a trimmed PDB (motif region +
substrate/partner + cofactor where available) and pre-styled `view.pml` / `view.cxc`
(grey protein, motif coloured as sticks with CPK heteroatoms + label, substrate/cofactor
shown, H-bond dashes). Open a script in PyMOL/ChimeraX, orient, save. Every structure
below was numbering-verified (catalytic residues / probe alignment) before mapping.
Manifests: `structures/<class>/*_motif_structures.csv`. The complete catalogue of
structures used:

**DNA — `structures/DNA/` (2 structures, both substrate + cofactor):**
- **`1MHT`** — *M.HhaI* (from *Haemophilus haemolyticus*), C5-cytosine DNA-MTase. X-ray,
  2.6 Å. Covalent ternary complex: protein + double-stranded **DNA** with the flipped
  5-fluorocytosine target base (`C36`) + **SAH** (the demethylated cofactor). Native
  numbering (Cys81 catalytic, Glu119 Motif VI). **Carries motifs 1–6.**
- **`1G38`** — *M.TaqI* (*Thermus aquaticus*), N6-adenine DNA-MTase. X-ray, 2.0 Å.
  protein + **DNA** + AdoMet-analogue cofactor (`NEA`) + N6-methyladenine product (`6MA`).
  Asn105 catalytic (NPPY). **Carries motifs 7–8.**

**RNA — `structures/RNA/` (3 PDB + 1 AlphaFold; mixed states):**
- **`2ZZM`** — *aTrm5* (*Methanocaldococcus jannaschii*), tRNA m1G37 MTase. X-ray.
  protein + cognate **tRNA-Leu** + **SAM**. *substrate + cofactor.* **Motifs 3, 6, 7.**
- **`9U48`** — *METTL16* (*S. pombe*), U6 snRNA m6A MTase. X-ray. catalytic domain + **RNA**
  + **SAM** (UniProt−13 offset). *substrate + cofactor.* **Motif 2.** (NB the related
  `9M86` is only the C-terminal KA1 domain and does **not** contain motif 2 — we used 9U48.)
- **`3M6U`** — *RsmF*/RsmB-family (*Thermus thermophilus*), 16S rRNA m5C MTase. X-ray.
  **apo** (no cofactor, no RNA) — motif-on-fold only. **Motifs 1, 4, 5.**
- **`AF-Q2GDL7`** — *RlmE* (*Ehrlichia sennetsu*), 23S rRNA 2′-O-MTase. **AlphaFold DB v6**
  predicted monomer, no ligands. **Motif 8.**

**Protein — `structures/Protein/` (2 complexes + 4 AlphaFold):**
- **`2NXN`** — *PrmA* (*Thermus thermophilus*), ribosomal-protein-L11 MTase. X-ray.
  enzyme (chain A) **+ its protein substrate L11** (chain B) — but **no cofactor** (apo of
  SAM). The only true enzyme-substrate complex in the set. **Motifs 1, 7.**
- **`4GQB`** — *PRMT5:MEP50* (human), protein-arginine MTase. X-ray. PRMT5 (A) + **MEP50
  partner/scaffold** (B) + peptide (C) + SAM-analogue cofactor (`0XU`). The "other protein"
  is a **partner, not the methylation substrate**. **Motifs 4, 6.**
- **`AF-Q6C4Y5`** (SET-domain histone-lysine MTase, *Yarrowia lipolytica*), **`AF-A7GSD9`**
  (malonyl-ACP O-MTase), **`AF-Q8DHV7`** and **`AF-Q9A9T7`** (HemK/PrmC release-factor
  glutamine MTases) — **AlphaFold v6** monomers, no ligands. **Motifs 2, 3, 5, 8.**

**Other — `structures/Other/` (all HOMOLOGY-mapped onto SAH-bound plant O-MTases):**
*The matched small-molecule MTases had zero SAM-bound structures, so each motif's PWM was
FIMO-matched onto a classic SAH-bound plant O-MTase co-crystal — a conserved feature shown
on a homolog, not the original protein. Figures must be captioned as such.*
- **`1KYZ`** — caffeic-acid 3/5-O-MTase (COMT) from alfalfa (*Medicago sativa*). X-ray.
  protein + **SAH** + **ferulic acid (`FER`) substrate**. *cofactor + substrate.*
  **Motifs 1, 3, 5, 6, 7, 8.**
- **`1M6E`** — salicylic-acid carboxyl MTase (SAMT, SABATH family) from *Clarkia breweri*.
  X-ray. protein + **SAH** + **salicylic acid (`SAL`)**. *cofactor + substrate.* **Motif 4.**
- **`3C3Y`** — PFOMT (phenylpropanoid/flavonoid O-MTase) from ice plant
  (*Mesembryanthemum crystallinum*). X-ray. protein + **SAH** (no substrate). *cofactor
  only.* **Motif 2.**

**IG-on-structure (Section 7), for reference — `structures/IG_<class>/`:**
`1YF3` (T4 Dam, X-ray, DNA + AdoHcy, *non-cognate*), `5WT1` (PaTrm5a, X-ray, tRNA + SAH,
*cognate* 4.4 Å), `5ZTZ` (NDUFAF7, X-ray, **apo**), `8PB3` (PsiM, X-ray, SAH + norbaeocystin
substrate, *cognate* 3.5 Å) — all four are experimental X-ray crystal structures.

---

## 9. Supporting analyses

- **ProstT5 self-attention** (`attribution/attention_analysis.py`, `attn_output/`): on the
  showcase proteins, attention does not globally track IG, but its *peak* attention lands
  on the catalytic/cofactor motifs (DPPY at 90th percentile, SAM loop at 80th) — the model
  internally flags the functional residues even though substrate discrimination is
  distributed.
- **Showcase per-residue IG** (`attribution/ig_output_logit/`): five reference enzymes,
  one per class, used to validate IG against known catalytic residues (e.g. C5-cytosine
  Cys160, Dam DPPY).
- **Deferred — ProstT5 SAE** (`interprot/PROSTT5_SAE_SCOPE.md`): scope for an InterProt-
  style sparse autoencoder to *name* the internal ProstT5 features; not yet run.

---

## Pipeline / file index

```
1_build_corpus.py / 2_make_reps.py / 3_domain_enrichment.py   stages 0–2
run_dedup.sh / run_domain_enrichment.sh / run_streme.sh        SLURM drivers
streme_<class>/                STREME motifs (txt, html, logos)
fimo_<class>/                  motif → rep matches (no-motif derivation)
domain_enrichment.csv          InterPro enrichment (Part A)
attribution/
  integrated_gradients.py      core IG (mean-pool-aware, logit space, completeness)
  no_motif_attribution.py      IG over all 3,023 no-motif reps + aggregate trends
  no_motif_attr_output/        per_sequence.csv, aa_enrichment, peakedness, top_windows
  attention_analysis.py        ProstT5 attention vs IG
  ig_output_logit/             showcase per-residue IG
  ig_q9v2g1/ ig_q54s83/ ig_p0dpa9/   per-protein IG for the structure-mapped cases
structures/
  <class>/                     motif-on-structure (PDB + view scripts) per class
  IG_<class>/                  IG-on-structure (B-factor PDB + view scripts) per class
interprot/                     deferred SAE scope + the InterProt paper
```

Reproduction note: large intermediates (`corpus.csv`, `corpus_all.fasta`, `clust_*`,
`tmp_mmseqs/`, `*_prostt5.h5` per-residue embeddings) are regenerable from the scripts and
are *not* shipped in the journal_article copy; structures, motif outputs, attribution
results, and figures are.
