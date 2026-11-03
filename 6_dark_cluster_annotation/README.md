# Stage 6 — Annotating the dark clusters: CLEAN + Horizyn cross-validation

Stage 5's global UMAP surfaces 33 tight, family-specific "Metabolite"-substrate clusters
(14 Actinomycetota, 19 plant; `dark33_members.csv`, 3,478 proteins total) that have
**zero reviewed/SwissProt members** — no experimentally characterized enzyme anchors
them, so the classifier's substrate-type call is the only annotation they have. This
stage asks, for each cluster: *is there any independent computational evidence for
what it actually does, and is it even a genuine methyltransferase?*

Two orthogonal, unsupervised methods are cross-validated against each other:

1. **CLEAN** (`clean/`) — EC-number assignment by nearest-neighbor distance (the
   "d-score") to per-EC embedding centroids learned by contrastive training on
   SwissProt. Run via `run_clean_dark33.sh` (maxsep decision rule) and
   `run_clean_dark33_pvalue.sh` (p-value rule) against `clean/dark33.fasta`. Lower
   d is more confident; **reliable / marginal / unreliable** tiers follow CLEAN's
   own published cutoffs. `clean/parse_clean_results.py` aggregates the raw
   per-protein calls (`dark33_maxsep.csv`) to a per-cluster dominant EC
   (`clean_dark33_per_cluster_ec.csv`) — this step is fully mechanical. Turning an
   EC number into a human-readable enzyme name and an MTase-yes/no judgment is
   **manual EC-nomenclature lookup**, captured as-curated in
   `clean/clean_dark33_cluster_summary.csv`.

2. **Horizyn** (`horizyn/`) — a pretrained dual-encoder contrastive model that
   embeds a *reaction* (RDKit+/DRFP fingerprints → 512-D) and a *protein*
   (ProtT5-XL → 512-D) into a shared space; cosine similarity ranks candidate
   reactions for a protein. `run_horizyn_predict_full.py` screens **all 3,478 dark
   proteins against all 36,552 Rhea reactions** (every EC class, not just
   SAM-methylation) — a deliberately hard, unbiased test. As a validation of the
   approach itself: with no hint that these are methyltransferases, Horizyn ranks
   a true SAM-methylation reaction as the single best match (rank 1 of 36,552) for
   **85.7%** of proteins (2,981 of 3,478), against a 4.7% random-chance baseline.
   Rhea lists every reaction twice, once per direction, and Horizyn scores the two
   entries almost identically (embedding cosine median 0.998), so a reaction counts
   as a SAM-methylation hit in either direction; counting the forward entry only
   gives 44.0% against 2.3% chance. The top-100 rate (98.6%) is not informative,
   because the chance level is already 91% to 99%. `horizyn/validation_by_direction.py`
   reproduces these numbers; per-protein ranks are in `horizyn/horizyn_full_validation.csv`.

   Because Rhea encodes one real chemistry as many near-duplicate reactions (e.g.
   one reaction per isoprenoid chain length for menaquinone), raw top-1 reaction
   IDs fragment a single activity across dozens of "distinct" hits. To fix this,
   `horizyn/cluster_rhea_reactions.py` clusters all 36,552 reactions in Horizyn's
   own 512-D embedding space with **UMAP (`n_neighbors=15, min_dist=0.1,
   metric=cosine`, `random_state=42`) → HDBSCAN (`min_cluster_size=5,
   min_samples=3`)** — the same embedding-clustering convention used everywhere
   else in this project (`../5_umap/global_reumap_prostt5.py`), with smaller
   HDBSCAN parameters because reaction "families" (chain-length variants of one
   chemistry) are much smaller than enzyme families. Output:
   `horizyn/rhea_reaction_clusters.csv` (per-reaction family label),
   `horizyn/reaction_family_summary.csv` (per-family size/orphan-rate/sample
   equations).

   `horizyn/remap_to_reaction_clusters.py` then re-labels each dark protein-cluster
   by its members' **dominant reaction family** (each protein's individual rank-1
   reaction, grouped into that reaction's family) rather than a single raw
   reaction. **HDBSCAN's noise label (`-1`) is not a real family** — it's the
   scattered leftover pile (15.9% of all reactions, but a protein-cluster's vote
   share going to `-1` should be read as "no coherent family," never as a
   describable chemistry). `horizyn/protein_clusters_family_final.csv` is the
   corrected final table: it determines each cluster's `coherent_dominant_family`
   from votes *excluding* `-1`, while separately reporting `pct_noise_votes` — an
   earlier version (not shipped here) mistakenly treated `-1` as a legitimate
   family and reported a misleading "98%+ orphan" for what was really just noise.

`family_descriptions.csv` gives a short, hand-verified (against real Rhea
equations pulled from the Rhea SPARQL endpoint) description of each reaction
family's chemistry, and flags families that are **not methylation at all**
(oxidoreductase / glycosyltransferase chemistry — see Headline results).

## Building the final table

`build_table2.py` joins `clean/clean_dark33_cluster_summary.csv` +
`horizyn/protein_clusters_family_final.csv` + `family_descriptions.csv` →
`table2_dark_clusters.csv`, reproducing the data behind the manuscript's Tables S1
(Actinomycetota) and S2 (Plant) — CLEAN's and Horizyn's independent calls side by
side for all 33 clusters.

## Headline results

- **Two strong, novel-methylation discovery leads**, where the methods agree
  something real and uncharacterized is happening without agreeing on the exact
  chemistry: **c758** (Horizyn: demethylmenaquinol methylation / vitamin K2
  biosynthesis, 61.7% of members individually vote for an orphan reaction in this
  family as their single best match of 36,552; CLEAN: thiol S-methyltransferase,
  marginal confidence) and **c181** (Horizyn: Co-precorrin methylation / cobalamin
  (vitamin B12) biosynthesis, 98.3% convergence, 0% noise; CLEAN: demethylspheroidene
  O-methyltransferase, unreliable confidence).
- **Four clusters both methods agree are NOT methyltransferases**, mis-swept into
  the dark SAM-MTase set: **c57, c200, c635** (100% unanimous vote for a
  spermidine/spermine synthase aminopropyl-transfer reaction; CLEAN independently
  calls all three "No" for MTase), and **c542** (70.6% vote for NADPH-dependent
  oxidoreductase chemistry; CLEAN's nearest neighbor is alpha-onocerin synthase, a
  terpene cyclase — different specific identity, same "not a methyltransferase"
  conclusion).

## Outputs

| File | Contents |
|---|---|
| `dark33_members.csv` | Ground-truth membership for all 33 dark clusters (Entry, Organism, cluster) — shared with Stage 5's dark-cluster grids and used throughout this stage. |
| `clean/clean_dark33_per_cluster_ec.csv` | Per-cluster dominant EC call, mechanically aggregated from CLEAN's raw per-protein output. |
| `clean/clean_dark33_cluster_summary.csv` | Curated per-cluster enzyme name / MTase yes-no / confidence tier (manual EC lookup). |
| `horizyn/protein_clusters_family_final.csv` | Per-cluster dominant reaction family, share, and noise fraction. |
| `horizyn/reaction_family_summary.csv` | Per-reaction-family size, orphan rate, sample equations (all 2,441 families incl. noise). |
| `horizyn/horizyn_full_validation.csv` | Per-protein rank/percentile of the best-scoring true SAM-methylation reaction (validation stat source). |
| `table2_dark_clusters.csv` | Final combined table — the data behind manuscript Tables S1 and S2. |

**Archived on Zenodo, not in git** (large): `horizyn_dark_scores.npz` (the full
3,478×36,552 score matrix, ~425 MB), `horizyn_dark_top50.csv` (~9 MB per-protein top-50) and
`horizyn_rhea_reaction_embeddings.npz` (the 512-D embedding of every Rhea reaction); `make data`
puts them in `data/`. Fig 8 and `horizyn/validation_by_direction.py` read them. To regenerate the scores, run
`horizyn/run_horizyn_predict_full.py` (requires the Horizyn `inf` checkpoint and a ProtT5-XL embedding of the
3,478 dark proteins via `horizyn/embed_prott5.py`).

## Caveats

- The reaction-family-clustering step groups reactions by embedding similarity,
  not by verified shared mechanism — a family's "chemistry" description is a
  human read of its sample equations, not a guarantee every member shares the
  exact same enzyme activity.
- CLEAN's EC-centroid distance is a nearest-neighbor call, not a functional
  proof; "reliable/marginal/unreliable" reflects distance-to-centroid, not
  experimental validation.

## Paths and environment variables

- Both tools run in their own conda environments (not the repository's `environment.yml`); the versions used are in
  `clean/environment_clean.yml` and `horizyn/environment_horizyn.yml`, and the drivers activate `clean` / `horizyn`
  (`CLEAN_ENV` / `HORIZYN_ENV` rename them).
- `clean/run_clean_dark33*.sh` need `CLEAN_APP_DIR` (the `app/` directory of a CLEAN checkout); the Horizyn
  scripts and `horizyn/run_*.sh` need `HORIZYN_DIR` (a Horizyn checkout holding `checkpoints/` and `configs/`). The Horizyn
  python scripts import the `horizyn` package from the folder they sit in, so copy `horizyn/*.py` into `$HORIZYN_DIR` before
  running the drivers (they `cd` there); outputs are written next to the scripts, so copy them back into this folder.
- `horizyn/inputs/` holds the Rhea reaction table (`rhea-reaction-smiles.tsv`) and the SAM-methylation reaction
  list the Horizyn scripts read; point `SMSTPRED_HORIZYN_INPUTS` elsewhere to override.
- `dark33.fasta` and `dark33_members.csv` are read from this folder (`SMSTPRED_STAGE6` overrides).

## Running it on a release (no manual steps)

`workflow/run_refresh.sh` runs this stage after the clusters of Stage 5, as five jobs (`--from dark_inputs --to annotate`):

| Job | Script | What it does |
|---|---|---|
| `dark_inputs` | `run_dark_inputs.sh` -> `make_dark_inputs.py` | the release's dark clusters (`<release>/umap/dark_clusters.csv`) -> `<release>/dark_annotation/<name>.fasta`, `<name>_members.csv`, `dark_name.txt` (`<name>` = `dark<number of clusters>`) |
| `clean` | `run_clean_release.sh` | CLEAN maxsep + p-value on that FASTA (`clean/`) |
| `horizyn` | `run_horizyn_release.sh` | Horizyn against all 36,552 Rhea reactions; the run folder `<release>/dark_annotation/horizyn_run/` links the checkout's pieces, so the checkout is never written to |
| `annotate` | `run_tables.sh` | CLEAN summary (tiers at d = 6 / 8, enzyme names from ExPASy ENZYME), Horizyn family table, `table_dark_clusters.csv`, Tables S1 and S2 (`tables/`), then Fig 8 and Fig 7 with Fig 8's clusters outlined |

CLEAN and Horizyn are external checkouts: their paths and conda environments are the `annotation` block of `refresh_config.json` (`CLEAN_APP_DIR`, `HORIZYN_DIR` override them).
`0_uniprot_pull/run_reference_tables.sh` (job `ref_tables`) downloads ExPASy ENZYME, `rhea2ec.tsv` and `rhea-directions.tsv` into `<release>/reference_tables/`.
What stays manual: the biological-context notes of the archived Table 2 ("-> anthracyclines") and the plain-language chemistry of reaction families beyond the 17 in `family_descriptions.csv` (new families carry a `[draft]` label made from their SAM-methylation flags).
Tests: `make annotation-test`.
