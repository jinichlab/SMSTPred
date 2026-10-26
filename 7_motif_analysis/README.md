# Sequence-Motif & Structural Analysis of SAM-MTase Substrate Specificity

This directory contains the **motif / interpretability** stage of the project: a set of
analyses that ask **what sequence and structural features distinguish methyltransferase
(MTase) substrate types** (DNA / RNA / Protein / Other-small-molecule). It combines two
complementary lines of evidence:

1. **Discriminative linear-motif discovery** (STREME, scanned with FIMO) — what short sequence
   motifs separate each substrate class, and how many enzymes carry none;
2. **Structural mapping** — placing the class-specific motifs onto experimental structures of
   real motif-carrying enzymes and measuring their distance to the cofactor and the substrate.

Together they show that the discriminative motifs recover the canonical catalytic and
cofactor signatures of each class, while a substantial fraction of enzymes (about 12 %)
carry no discriminative motif at all and sit scattered inside motif-bearing families — the
substrate signal is distributed rather than a short linear motif.

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

### The class-specific motifs — mapped structure and structural placement

Of the 32 discovered motifs, **26 are class-specific** (own-class carrier rate ≥ 2× the
pooled-background rate in the FIMO cross-scan over all 118,816 proteins; the other 6 are
pan-methyltransferase SAM-binding / Rossmann-core motifs — see §2.1). For each of the 26 we
select, **from FIMO evidence rather than by hand**, the best in-class carrier whose
experimental structure has **both** a SAM/SAH cofactor and the substrate bound, locate the
motif by its FIMO `matched_sequence` (author numbering), and measure the minimum
heavy-atom distance from the motif to the cofactor and to the substrate. A motif is **near**
a ligand at **≤ 4 Å** (direct van der Waals / H-bond contact); placement is **near substrate**,
**near SAM**, **near substrate/SAM** (active site), or **distal** (> 4 Å from *both*).
A geometry gate rejects structures where the ligand sits on a different
chain/subunit (e.g. a ribosome-assembly complex), falling back to the next carrier.

**#** is a sequential index within class; **Motif** is the original STREME ID (gaps are the
6 non-specific motifs removed above). Per-motif PDB, chain, range and distances are in
`motif_placements.csv`; the folder for each is `structures/<Class>/<class>_motif<STREME#>_<PDB>/`.
See `structures/rebuild_from_carriers.py` for the pipeline.

**DNA**

| # | Motif | Consensus | Structure | Placement |
|---|---|---|---|---|
| 1 | DNA-1 | `IGGPPCQGFSRAG` | M.HhaI `10MH` — X-ray; **SAH + DNA** | **near substrate/SAM** — SAM 3.7, sub 2.6 Å |
| 2 | DNA-2 | `PFAGSGGT` | CMT3 chromomethylase `7UBU` — X-ray; **SAH + DNA** | **near SAM** — SAM 2.9, sub 6.0 Å |
| 3 | DNA-3 | `YGVPQRRERVF` | M.HhaI `10MH` — X-ray; **SAH + DNA** | **near substrate** — SAM 8.3, sub 2.7 Å |
| 4 | DNA-4 | `LTPREAARLQGFPDDYTF` | CMT3 chromomethylase `7UBU` — X-ray; **SAH + DNA** | **distal** — SAM 5.1, sub 5.8 Å |
| 5 | DNA-5 | `KQIGNAVPPPVARAJA` | CMT3 chromomethylase `7UBU` — X-ray; **SAH + DNA** | **near substrate/SAM** — SAM 2.7, sub 3.5 Å |
| 6 | DNA-6 | `FLLENVPGL` | M.HhaI `10MH` — X-ray; **SAH + DNA** | **near substrate** — SAM 5.0, sub 3.0 Å |
| 7 | DNA-7 | `YTBPPYNT` | M.TaqI `1G38` — X-ray; **NEA + DNA** | **near substrate/SAM** — SAM 2.9, sub 2.8 Å |
| 8 | DNA-8 | `GZFFTPREI` | EcoKI (Type I R-M) `2Y7H` — cryo-EM; **SAM + DNA** | **near substrate/SAM** — SAM 2.5, sub 2.7 Å |

**RNA**

| # | Motif | Consensus | Structure | Placement |
|---|---|---|---|---|
| 1 | RNA-1 | `DRVLLDPPCSGLGT` | NSUN6 `5WWR` — X-ray; **SFG + RNA** | **near substrate/SAM** — SAM 2.4, sub 2.4 Å |
| 2 | RNA-2 | `RJVYSTCS` | NSUN6 `5WWR` — X-ray; **SFG + RNA** | **near substrate** — SAM 5.9, sub 2.6 Å |
| 3 | RNA-4 | `LQDASSMLPALALDPQP` | NSUN6 `5WWR` — X-ray; **SFG + RNA** | **near substrate** — SAM 5.3, sub 2.7 Å |
| 4 | RNA-5 | `TGGFFVAKL` | NSUN6 `5WWR` — X-ray; **SFG + RNA** | **near substrate** — SAM 6.6, sub 3.7 Å |
| 5 | RNA-6 | `ENEAVVAAALK` | NSUN6 `5WWR` — X-ray; **SFG + RNA** | **distal** — SAM 10.8, sub 4.1 Å |
| 6 | RNA-7 | `KKKKKKRK` | aTrm5 `2ZZM` — X-ray; **SAM + RNA** | **distal** — SAM 9.1, sub 5.2 Å |
| 7 | RNA-8 | `RSRAAFKLJZL` | TrmJ/FTSJ1 `8Y2O` — cryo-EM; **SAH + RNA** | **near substrate/SAM** — SAM 3.3, sub 2.8 Å |

**Protein**

| # | Motif | Consensus | Structure | Placement |
|---|---|---|---|---|
| 1 | Protein-2 | `GGTGLKVWEAAJVLA` | METTL21D/VCPKMT `7OAT` — X-ray; **SAH + protein substrate** | **near substrate/SAM** — SAM 3.3, sub 2.8 Å |
| 2 | Protein-3 | `AADVVYDP` | METTL21D/VCPKMT `7OAT` — X-ray; **SAH + protein substrate** | **near substrate/SAM** — SAM 2.5, sub 2.6 Å |
| 3 | Protein-4 | `VDIJVSEWMGYFLLNES` | PRMT1 `1OR8` — X-ray; **SAH + protein substrate** | **near substrate/SAM** — SAM 3.5, sub 2.8 Å |
| 4 | Protein-5 | `THWKQTVFYL` | PRMT1 `1OR8` — X-ray; **SAH + protein substrate** | **near substrate** — SAM 7.8, sub 2.8 Å |
| 5 | Protein-6 | `VLIPRPETE` | PrmC/HemK `2B3T` — X-ray; **SAH + protein substrate** | **near substrate/SAM** — SAM 2.8, sub 2.6 Å |
| 6 | Protein-7 | `AFGTGTHPTTRLCLEAL` | PrmA (L11 MTase) `3CJQ` — X-ray; **SAH + protein substrate** | **near substrate/SAM** — SAM 2.6, sub 2.6 Å |
| 7 | Protein-8 | `SNPPYIPT` | PrmC/HemK `2B3T` — X-ray; **SAH + protein substrate** | **near substrate/SAM** — SAM 3.4, sub 2.7 Å |

**Other**

| # | Motif | Consensus | Structure | Placement |
|---|---|---|---|---|
| 1 | Other-5 | `GGKERTEEEWRAL` | norcoclaurine 6-OMT `5ICE` — X-ray; **SAH + small-molecule substrate (`2H4`)** | **distal** — SAM 10.3, sub 6.5 Å |
| 2 | Other-6 | `LHDWPDEDCVKIL` | GedA OMT `7WH9` — X-ray; **SAH + small-molecule substrate (`EMO`)** | **near substrate/SAM** — SAM 3.3, sub 2.8 Å |
| 3 | Other-7 | `HDFFTPQPVKGADVYL` | MpaG' OMT `8XTE` — X-ray; **SAH + small-molecule substrate (`A1LWB`)** | **near SAM** — SAM 3.0, sub 6.9 Å |
| 4 | Other-8 | `GIVQDLPEVIA` | scoulerine 9-OMT `6I6K` — X-ray; **SAH + small-molecule substrate (`SLX`)** | **near SAM** — SAM 2.6, sub 9.2 Å |

**Summary.** All **26** class-specific motifs map onto an **experimental** structure of a
real MTase that carries the motif (X-ray; cryo-EM for the FTSJ1 panel), with **both** SAM/SAH
and substrate bound and the placement measured directly — no homology models, no AlphaFold-only
mappings, no hand-picked ranges. The motifs are overwhelmingly **active-site elements**:
13 sit at the substrate–cofactor junction (both ligands ≤ 4 Å), 6 contact only the substrate, 3 only the SAM pocket,
and **4** are **distal** — > 4 Å from *both* ligands: RNA-6 (`ENEAVVAAALK`), RNA-7 (`KKKKKKRK`), DNA-4 (`LTPREAARLQGFPDDYTF`), Other-5 (`GGKERTEEEWRAL`).
Of these, Other-5, RNA-7 and RNA-6 sit well off the cofactor (> 9 Å from SAM) on the
substrate-recognition face, while DNA-4 lies ~5–6 Å from both in the active-site cleft
(junction-adjacent, in second shell rather than direct contact). The 6 non-specific
motifs (RNA-3, Protein-1, Other-1/2/3/4) are excluded as pan-MTase SAM-binding / Rossmann-core
elements shared across classes.


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
signature).

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

### Where the no-motif proteins sit on the embedding landscape

`no_motif_umap.py` / `no_motif_umap_interactive.py` project motif coverage onto the
Stage-5 UMAP (coords reused, no recompute). Two-level labeling: each 50%-id cluster's
motif status (its rep, from FIMO) is propagated to **all 106,608 members** via
`clust_cluster.tsv`, then judged per UMAP (HDBSCAN) cluster.

- **A — near-motif-dark families:** only **12** UMAP clusters (about 960 proteins) are
  at least 90 % motif-less (91–99 %), mostly **RNA**; **no cluster is entirely motif-free**.
- **B — all no-motif proteins:** the **12,866** motif-less proteins (12.1 % of all 106,608
  class members) are **mostly scattered inside motif-bearing families**, i.e. divergent
  members rather than novel groups. (Consistent with the out-of-domain / fragment enrichment
  above.)

Outputs: `no_motif_umap.png` (A|B panels), `no_motif_umap_clusters.csv` (the 12 clusters),
and `no_motif_umap_interactive.html` — a `datamapplot` view with an **A/B toggle**, hover,
search, and click-to-UniProt.

---

## 6. Structure selection — from FIMO carriers, not by hand

Every structure panel in the table above is selected automatically from FIMO evidence
(`structures/rebuild_from_carriers.py`). For each of the 26 class-specific motifs we take the
best in-class carrier whose experimental structure has **both** a SAM/SAH cofactor and the
substrate bound (ligand characterization by `scan_pdb_ligands.py` over the cross-scan), locate
the motif by its FIMO `matched_sequence` in author numbering, measure motif→cofactor and
motif→substrate distances, and classify the placement. A geometry gate rejects structures where
the ligand sits on a different chain/subunit (e.g. a ribosome-assembly complex) and falls back
to the next carrier. Per-motif PDB, chain, residue range, matched sequence, distances and
placement are in `structures/motif_placements.csv`; each subset PDB + `view.pml` + `view.cxc`
is in `structures/<Class>/<class>_motif<N>_<PDB>/`.

**Rendering.** Every `view.pml` / `view.cxc` draws the protein grey (`#b0b0b0`), the motif in
its class colour with CPK side chains, any protein/nucleic substrate in wheat/tan, a
small-molecule substrate in **gold**, and the **SAM/SAH cofactor (and analogs `0XU`/`NEA`/`QVR`/
`SFG`) with warm-brown carbons (`#b5651d`)** so the cofactor never reads as protein grey.

*(Supersedes the earlier hand-curated per-PDB catalogue and the §10/§10.1 upgrade searches.)*

---

## 7. Structure due-diligence (superseded)

Earlier versions selected one representative structure per motif and then ran two manual
"upgrade" searches — a cluster-member search and a FIMO-carrier search — to move
motifs off apo/AlphaFold structures onto substrate/cofactor-bound ones. Both are **superseded**
by the systematic pipeline in §6: the FIMO cross-scan (`fimo_cross/`, all 32 motifs × 118,816
proteins), `scan_pdb_ligands.py`, and `rebuild_from_carriers.py` now pick the best
substrate+SAM-bound in-class carrier for every motif directly and re-measure placement. This
also corrected several mis-mappings from the manual era — e.g. RNA-2 had been placed on METTL16
`9U48`, which does not carry the motif at all; it is now on NSUN6 `5WWR` at its catalytic loop.

---

## Pipeline / file index

```
1_build_corpus.py / 2_make_reps.py / 3_domain_enrichment.py   stages 0–2
run_dedup.sh / run_domain_enrichment.sh / run_streme.sh        SLURM drivers
streme_<class>/                STREME motifs (txt, html, logos)
fimo_<class>/                  motif → rep matches (no-motif derivation)
domain_enrichment.csv          InterPro enrichment (Part A)
structures/
  <class>/                     motif-on-structure: each motif folder = full-structure PDB
                               + view.cxc + view.pml (grey fold, motif highlighted); built by
  make_context_views.py          self-contained context views (27 motifs, from source PDBs)
  make_context_views_upgrades.py the 5 bound-structure upgrades (§7, from _upgrade_src/)
find_bound_structures.py       §7 cluster-member → UniProt→PDB → bound-structure search
verify_upgrade_motifs.py       FIMO-verify each upgrade structure carries its motif
no_motif_prep.py / no_motif_manifest.csv   the 3,023 no-motif family reps (Entry, class, length)
no_motif_umap.py / _interactive.py   motif coverage on the embedding UMAP (§5)
interprot/                     scope note for planned (not yet run) per-residue attribution and SAE work
```

Reproduction note: large intermediates (`corpus.csv`, `corpus_all.fasta`, `clust_*`,
`tmp_mmseqs/`, `*_prostt5.h5` per-residue embeddings) are regenerable from the scripts and
are *not* shipped in the journal_article copy; structures, motif outputs, and figures are.
