# Stage 6 — Annotating the dark clusters: CLEAN + Horizyn cross-validation

Stage 5's global UMAP surfaces 33 tight, family-specific "Other"-substrate clusters
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
   **44%** of proteins, and within the top 100 for **98.6%**, against a 2.3%
   random-chance baseline (median rank 2) — see `horizyn/horizyn_full_validation.csv`.

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
`table2_dark_clusters.csv`, reproducing the data behind the manuscript's Table 2
(a: Actinomycetota, b: Plant) — CLEAN's and Horizyn's independent calls side by
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
| `dark33_members.csv` | Ground-truth membership for all 33 dark clusters (Entry, Organism, cluster) — shared with Stage 5's Fig. 6C/D and used throughout this stage. |
| `clean/clean_dark33_per_cluster_ec.csv` | Per-cluster dominant EC call, mechanically aggregated from CLEAN's raw per-protein output. |
| `clean/clean_dark33_cluster_summary.csv` | Curated per-cluster enzyme name / MTase yes-no / confidence tier (manual EC lookup). |
| `horizyn/protein_clusters_family_final.csv` | Per-cluster dominant reaction family, share, and noise fraction. |
| `horizyn/reaction_family_summary.csv` | Per-reaction-family size, orphan rate, sample equations (all 2,441 families incl. noise). |
| `horizyn/horizyn_full_validation.csv` | Per-protein rank/percentile of the best-scoring true SAM-methylation reaction (validation stat source). |
| `table2_dark_clusters.csv` | Final combined table — the data behind manuscript Table 2. |

**Not shipped** (large, regenerable): `horizyn_full_scores.npz` (the full
3,478×36,552 score matrix, ~425 MB) and `horizyn_full_topk.csv` (~9 MB
per-protein top-50) — regenerate with `horizyn/run_horizyn_predict_full.py`
(requires the Horizyn `inf` checkpoint and a ProtT5-XL embedding of the 3,478
dark proteins via `horizyn/embed_prott5.py`).

## Caveats

- The reaction-family-clustering step groups reactions by embedding similarity,
  not by verified shared mechanism — a family's "chemistry" description is a
  human read of its sample equations, not a guarantee every member shares the
  exact same enzyme activity.
- CLEAN's EC-centroid distance is a nearest-neighbor call, not a functional
  proof; "reliable/marginal/unreliable" reflects distance-to-centroid, not
  experimental validation.
