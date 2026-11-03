# Stage 0 — Refresh the UniProt inputs

Optional and upstream of everything else. The stages 1–7 run on the **archived** data (Zenodo). This stage pulls the **current**
UniProt data through the UniProt REST API so the analysis can be re-run on a newer release:

| Part | UniProt query (sent to `rest.uniprot.org`) | What it is |
|---|---|---|
| `sam` | `(reviewed:true) AND (xref:interpro-IPR029063)` | Swiss-Prot proteins of the SAM-dependent methyltransferase superfamily (Class I): the SAM training set before labelling |
| `up` | `(reviewed:true) AND (cc_catalytic_activity:*)` | Swiss-Prot enzymes with a catalytic-activity annotation: the broad UP comparison set before labelling |
| `proteomes` | `(proteome:UPxxxxxxxxx)` for every Proteome Id of [`data/selected_1259_organisms.csv`](../data/selected_1259_organisms.csv) | the 1,259 reference proteomes of the tree-of-life analysis (all their UniProtKB proteins) |

The queries reproduce how the archived sets are defined; the archived tables were checked for it: the 7,405 SAM proteins all carry
IPR029063 and are reviewed; the 201,339 UP proteins are all reviewed, all have catalytic-activity text (92 % also a Rhea ID) and none is
longer than 1,000 aa. The archived tables were additionally **labelled** (Metabolite / Protein / RNA / DNA) from the reaction annotation:
the second half of this stage (`bin_substrates.py`, below) does that deterministically from the pulled text, so the pulled sets are supersets
of the labelled ones (entries whose reactions cannot be typed stay `label_source = none`).

## Length: nothing is dropped at pull time

The archived analysis embedded and predicted **<= 1,000 residues** (the archived labelled training sets stop at 1,000 aa). For the refresh the
deployed ProstT5 model is trained and applied on **<= 2,000 aa** (decision of 2026-11-02): it needs only the sequence, and a 2,000 aa limit keeps
about 93 % of the long DNA methyltransferases of the 1,259 proteomes (DNMT1-, DOT1L- and CMT-type enzymes) that a 1,000 aa limit loses. The pull
keeps every length and records both views (`--length-view`, default 2,000; 0 = no split; `resplit` changes it later without downloading):

* `sam_mtases.tsv.gz`, `up_enzymes.tsv.gz`: all lengths. `*_le2000.fasta`: the proteins that get embedded.
* `reference_proteomes/<organism>.csv`: proteins <= 2,000 aa, **same file names and columns** as the archived tables
  (`Entry, Reviewed, Protein names, Organism, Length, Sequence, InterPro, Pfam`), so `4_phylo_trends/reference_tree.py` and
  `7_motif_analysis/1_build_corpus.py` read them unchanged (point `SMSTPRED_DATA` at a folder that has them as `reference_proteomes/`).
  `reference_proteomes_long/<organism>.csv`: the proteins longer than that (kept, not used downstream; mostly polyketide / non-ribosomal
  synthases and viral polyproteins, which would need domain-level embeddings).
* `proteome_summary.csv`: per proteome the counts of both views (all, <= / > the cutoff, reviewed, IPR029063).

Methods that cannot take proteins over 1,000 aa are scored on the <= 1,000 aa subset (filter on the `Length` column): SaProt (1,024-token limit
in `embed_saprot.py`), and ESM-IF and the Foldseek baseline, which need a structure (AlphaFold DB has no full-length model above about 2,700 aa
and lacks some below it; of the 462 reviewed IPR029063 proteins over 1,000 aa, 196 have one).

Why the split is explicit: in the archived `reference_proteomes/` the 146 organisms re-pulled on 2026-09-08 (to add the `Reviewed` column) still
contain their proteins longer than 1,000 aa (24,408 proteins, 3.3 % of those files), while the other 1,113 files do not.
`reference_tree.py` counts every row of a file as "proteome size", so the proteome sizes and SAM counts of those 146 organisms (Fig 5) are
about 3 % larger than those of the others. The refreshed tables apply one cutoff to every organism. (Proteome sizes and SAM counts at a
2,000 aa cutoff are therefore not comparable with the archived Fig 5 table, which is <= 1,000 aa for 1,113 of the organisms.)

## Output and provenance

Everything goes to `$UNIPROT_PULL_DIR/release_<UniProt release>/` (default `<repo>/data/uniprot_pull/`, ignored by git); the archived data and
the committed tables are never touched. The UniProt release is read from the API (`X-UniProt-Release`), put in the folder name, and
checked on every response: if UniProt switches release during a pull the pull stops (`--expect-release 2026_03` makes the intended
release explicit). `pull_manifest.json` records the release and date, the queries and fields, row counts, the SHA-256 of every output and
the repository commit; a copy is kept in [`pull_manifests/`](pull_manifests) so it can be committed.

## Run

```bash
make pull                                   # = bash 0_uniprot_pull/run_pull.sh; on SLURM: sbatch 0_uniprot_pull/run_pull.sh
make pull-smoke                             # 3 proteomes + a few hundred rows, to check the API is reachable
make pull-compare                           # what changed against the archived data (also run at the end of run_pull.sh)

UNIPROT_PULL_DIR=/path/to/refresh PULL_ARGS="--expect-release 2026_03" sbatch run_pull.sh
PULL_WHAT=proteomes sbatch run_pull.sh      # one part only (sam | up | proteomes | all); the proteome part resumes after an interruption
PULL_WHAT=resplit PULL_ARGS="--length-view 1500" sbatch run_pull.sh    # another length cutoff for the newest pull, no download
```

| File | What it does |
|---|---|
| `pull_uniprot.py` | the three pulls (`sam`, `up`, `proteomes`, `all`), manifest, release pinning |
| `uniprot_api.py` | the REST client: rate limit, retries with backoff, release tracking, streamed / paged TSV downloads with row-count checks |
| `compare_to_previous.py` | new vs archived: SAM and UP proteins added / removed / with a changed sequence; per proteome the IPR029063 proteins vs the 112,019 predicted ones; Fig 5 proteome sizes and SAM counts |
| `run_pull.sh` | SLURM / bash driver |
| `substrate_binning.py`, `bin_substrates.py`, `test_substrate_binning.py`, `run_bin.sh` | substrate typing: rules, table driver, golden cases, SLURM / bash driver (section below) |
| `fetch_rhea_compounds.py`, `rhea_compounds.tsv.gz` | Rhea's family of every compound (polypeptide / polynucleotide / polymer / small molecule), the table the typing looks names up in; `polypeptide_tokens_reviewed.txt`: the polymer names reviewed so far |
| `pull_manifests/` | small provenance record of each pull (tracked) |

Network use: about 1,260 proteome requests plus two large downloads, spaced at least 0.2 s apart (4 parallel proteome workers).

## Substrate typing ("binning"): Metabolite / Protein / RNA / DNA

`bin_substrates.py` turns the pulled "Catalytic activity" text into the label tables the later stages read. The rules are in `substrate_binning.py`
(pure functions; hand-checked golden cases in `test_substrate_binning.py`; `make bin-test`). They replace the notebooks of the original labelling,
which decided "Protein" by matching words inside the `[...]` of a reaction (a new polymer name in a later annotation would have been missed).
The substrate is now typed from **Rhea's own classification of the compound**, not from its name: `fetch_rhea_compounds.py` downloads one table
(15,203 compounds, one SPARQL query, 160 kB, kept in this folder as `rhea_compounds.tsv.gz`) and the UniProt reaction text, which uses Rhea's
compound names verbatim, is looked up in it.

Per reaction block of the text:

1. **Which reactions.** `sam` set: reactions that have S-adenosyl-L-methionine as a reactant (methyl transfer, but also carboxy-SAM synthesis (CmoA)
   and aminopropyl transfer (nicotianamine / wybutosine synthesis); the kind is kept in the `transfer` column). `up` set: every reaction.
2. **The acceptor.** The reactants that are not SAM / water / protons / nucleotides / redox cofactors / **redox-partner proteins** (ferredoxin, thioredoxin,
   flavodoxin, plastocyanin, cytochromes, NADPH-hemoprotein reductase ...).
3. **Its class, from Rhea's family of the compound** (exact name lookup; `Na(+)(out)` is looked up as `Na(+)`):
   `SmallMolecule` -> Metabolite; `Polypeptide` (a protein-bound residue such as `L-histidyl-[protein]`, 1,268 compounds, whatever the polymer is called) -> **Protein**;
   `Polynucleotide` (712 compounds; every name carries an RNA or DNA word) -> RNA or DNA by that word; `Polymer` (264; Rhea mixes `RNA(n)`, `DNA(n)`,
   polyglutamylated proteins, sugar chains, polyphosphate and polyprenyls in it) -> the word rules. A name that is not in the table (UniProt annotated the reaction
   without a Rhea id: ubiquitin E2, protein phosphotransfer, plant O-methyltransferases) also goes to the word rules, which are: an RNA word (`tRNA`, `rRNA`, `RNA(n)`,
   token-bounded, so "external" is not RNA), a DNA word (`DNA`, `ssDNA`, `(deoxyribonucleotide)n`), a polymer bracket `residue-[name]` that begins with a letter and is not a sugar chain,
   polyphosphate or redox partner, or the word protein / peptide; otherwise Metabolite.
   **Rhea's family is not trusted blindly.** Its generic polypeptide class is its default for a compound with an R group, so 6 of the 1,268 polypeptide compounds are really
   nucleic-acid residues (`...selenouridine(34) in tRNA`, `an N(7)-methylguanosine in RNA`, `5-aminomethyl-dUMP in DNA`; 108 proteins in the 2026_03 pull) and two are small molecules
   (`a Fe(III)-siderophore`, `a Fe(II)-siderophore`): 8 mis-filings in all. A nucleic-acid word **outside** a bracket overrules the polypeptide family (inside a bracket, `[RNA ligase]` or `[DNA-directed RNA polymerase]`, it only names a protein);
   every such conflict is a line of `binning_rules_used_*.csv`; the two siderophores are `FAMILY_OVERRIDES` in `substrate_binning.py`.
   If Rhea corrects a filing, nothing has to be changed: the compound then arrives as Polynucleotide / SmallMolecule and is typed the same (checked on the real table: 0 class changes), the conflict rule simply stops
   firing, and each run prints which safeguards are still needed (`Rhea table check` in the report). Only a mis-filing of a *small molecule* (the siderophores) needs a manual `FAMILY_OVERRIDES` entry,
   and the "NEW polymer names" list is what shows it.
4. **Per reaction** a macromolecule replaces Metabolite (an aminoacyl-tRNA synthetase acts on tRNA, not on its amino acid); **per enzyme** the classes of its
   reactions are joined (14 multi-label enzymes in the archived SAM set).
5. **Evidence, gated per reaction.** A reaction is *curated* when its own evidence matches HAMAP / PROSITE-ProRule / RuleBase (a curator-written rule) or PubMed / Ref.N
   (a literature reference), the same patterns as the archived labelling; by-similarity (`ECO:0000250`) and automatic annotation do not count. The labels of an enzyme come **only
   from its curated reactions** (`label_source = reaction_curated`; the reactions that were left out are counted in `n_uncurated_ignored`). An enzyme none of whose SAM reactions is
   curated is typed from all of them but marked `reaction_uncurated`, and stays out of `in_benchmark` / `in_deployment` unless you add it with `--*-sources`.
   (The archived labelling gated per enzyme: one curated block let every block in. In the SAM set that admitted 88 uncurated SAM reactions and 40 enzymes whose SAM reactions are all uncurated.)
   Most "curated" evidence is rule-propagated, not experimental: of the SAM reactions behind the labels 86 % rest on a HAMAP / ProRule rule and 13 % on a literature reference;
   requiring experimental evidence would leave 691 of 7,553 enzymes.

What still needs a human, and how a refresh shows it: the **redox partners** (ferredoxin, thioredoxin, plastocyanin, flavodoxin, cytochromes, DsbA ...) are the one name list left, because
Rhea files them as ordinary polypeptides. A partner Rhea adds later would be typed Protein, so every run compares the polypeptide polymers that typed a Protein with
`polypeptide_tokens_reviewed.txt` (163 names, read through on 2026-11-02) and prints the ones that are new; a redox partner among them goes into `REDOX_PARTNER`,
`--accept-new-tokens` adds the rest to the list. `binning_text_fallback_names_*.csv` lists every name the word rules had to type (1.1 % of the SAM-set typings, 8.4 % of the UP-set).
6. **Weaker tier, SAM set only.** Enzymes without a typable reaction can be typed from their names (`label_source = name_rule`; RNA / DNA / Protein only, never
   Metabolite: "metabolite" is the absence of a macromolecule and a name cannot show an absence). It is scored against the reaction-typed enzymes in the report.

| Output (in the release folder) | Content |
|---|---|
| `labels_SAM.csv`, `labels_UP.csv` | the four 0/1 label columns (same names as the archived `labels_SAM_full.csv` / `labels_UP.csv`), `label_source`, `transfer`, `Protein_via_carrier` (Protein only because a carrier protein holds the substrate: ACP, peptidyl / sulfur / lipoyl carrier, PKS / NRPS / FAS domains; 70 SAM-set and 5,643 UP-set enzymes, the archived labelling counted them as Protein), `n_reactions` (the curated reactions that were used), `n_uncurated_ignored`, `RheaID`, `Acceptors` (every acceptor with the rule that typed it: `rhea:...` structural, `text:...` word rules), `in_benchmark`, `in_deployment` |
| `binning_rules_used_{sam,up}.csv`, `binning_text_fallback_names_{sam,up}.csv` | every rule that fired and how many proteins it typed (the review list); the names the word rules had to type |
| `binning_vs_archived_{sam,up}.csv` | the proteins whose label set differs from the archived labelling, with both sets and the acceptors |
| `binning_report.txt`, `binning_manifest.json`, `rhea_compounds.tsv.gz` | counts, structural share, new polymer names, agreement with the archive, accuracy of the name tier; input and Rhea-table sha256 and parameters; the Rhea table used |

`in_benchmark` = `Length <= --benchmark-max-len` (1,000) and `label_source` in `--benchmark-sources`; `in_deployment` = `Length <= --deployment-max-len` (2,000) and
`label_source` in `--deployment-sources`. Both source lists default to `reaction_curated`, which is what the archived labelling used.

### Check against the archived labelling (UniProt release 2026_03)

| | proteins in both | identical label sets | what the differences are |
|---|---|---|---|
| SAM set (`labels_SAM_full.csv`) | 7,404 of 7,405 (the other has no typable reaction any more) | **7,388 (99.78 %)** | 6 arsenite methyltransferases: archived Protein only because `[thioredoxin]` is in the reaction (the acceptor is arsenic triglutathione; now Metabolite); 6 (3 fibrillarins, 3 DNMT3A) whose archived extra Protein label came from a by-similarity reaction that the per-reaction gate ignores; 4 release changes (CARNMT1 and others) or the cypemycin N-terminal MTase, whose acceptor is a ribosomal precursor peptide (now Protein, a judgement) |
| UP set (`labels_UP.csv`) | 200,710 of 201,339 | **190,901 (95.11 %)** | about 4,970 archived Protein only because of a redox-partner protein; 1,246 DNA-directed RNA polymerase (archived DNA + RNA, the reaction text has RNA only); 459 autocleaving repressors (LexA, cI: archived DNA, the reaction is protein hydrolysis); 830 archived RNA from tRNA-dependent enzymes whose reaction text has no tRNA; the other ~2,300 are mostly an added Metabolite or Protein reaction and were not audited one by one |

Net effect on the UP classes (the 200,710 proteins in both): Protein 30,899 -> 25,332, DNA 5,689 -> 4,131, RNA 36,335 -> 35,505, Metabolite 133,243 -> 138,899. The differences above were counted before the Rhea-family lookup and the percentages moved by under 0.1 points when it replaced the word rules: the two routes agree.

Name tier on the SAM set, scored against the reaction labels: on the 7,722 curated enzymes it fires on 80 % with precision 99.9 % (RNA), 96.9 % (DNA), 99.8 % (Protein);
on the 850 by-similarity enzymes (less studied, the fairer estimate) it fires on 65 % with 91.7 / 100 / 99.5 %. A "methyltransferase with no macromolecule word -> Metabolite" rule
scores 100 % on the curated enzymes but is **not** used: the unlabelled remainder holds TFB1M (an rRNA dimethylase whose name hides it), SAMTOR (a SAM sensor, not an
enzyme) and NDUFAF5 (an arginine hydroxylase). The remaining 1,241 untyped SAM-set proteins are mostly not SAM methyltransferases (spermidine synthases, anamorsin, fatty acid synthase).

### Run

```bash
make bin-test                               # the golden cases, offline, seconds
make bin                                    # = bash 0_uniprot_pull/run_bin.sh; on SLURM: sbatch 0_uniprot_pull/run_bin.sh (about 2 minutes, CPU)
BIN_WHAT=sam sbatch run_bin.sh              # one set only (sam | up | all)
BIN_ARGS="--deployment-max-len 1500" sbatch run_bin.sh
BIN_ARGS="--benchmark-sources reaction_curated reaction_uncurated" sbatch run_bin.sh   # also count the by-similarity enzymes
BIN_ARGS="--refresh-rhea" sbatch run_bin.sh                 # download a new Rhea compound table first (network); run after every UniProt refresh
BIN_ARGS="--accept-new-tokens" sbatch run_bin.sh            # after reading the "NEW polymer names" list of the report
```

## Not part of this stage (yet)

Re-embedding the new sequences (`embeddings/`), re-clustering and re-training, and re-predicting the proteomes. The pulled tables carry the fields those steps
need (catalytic-activity text, Rhea IDs, EC, InterPro, Pfam, the entry and sequence dates and versions). `make verify` does not run this stage (it is a tier-C
step in `workflow/steps.py`).
