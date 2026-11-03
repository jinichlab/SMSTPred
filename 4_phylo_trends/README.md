# Stage 4 — Phylogenetic & taxonomic trends in predicted substrate usage

Turns the tree-of-life predictions (Stage 3) into organism-level and lineage-level trends:
how predicted MTase substrate usage varies across the tree and with host lifestyle.

## Scripts

| Script | Purpose |
|---|---|
| `reference_tree.py` | Builds a taxonomic **Newick tree** over the ~1,259 reference organisms and annotates each leaf with its per-substrate call composition for **iTOL** rendering (ring/bar tracks of DNA/RNA/Protein/Metabolite usage). |
| `protein_cellularity_violin.py` | Per-organism **Protein-substrate** MTase load (count or fraction) split across eukaryotic **cellularity groups** (e.g. metazoa / unicellular / plant), on the ProstT5 deployment predictions — testing whether protein-methylation capacity tracks multicellular complexity. |
| `run_reference_tree.sh` | SLURM driver. |
| `run_tree.sh` | **Refresh release, no iTOL:** `3_predictions/combine_predictions.py` -> `reference_tree.py --release` -> `figures/Fig3_tree_landscape/draw_tree.py`; writes `<release>/tree/` (iTOL files, `tree_stats.csv`, `fig3_tree.png/.svg/.pdf`). `make tree`. |
| `test_tree.py` | tests of the release-mode steps on synthetic data. |

In release mode every call of the ProtT5 ensemble counts (in the applicability domain or not); a ring value is the share of the organism's *predicted* SAM-MTases (<= 2,000 aa)
called a substrate type (a protein with two calls counts in both), while the bar and "% SAMs in proteome" use every IPR029063 protein. The tree is the NCBI taxonomy with unit branch lengths
(a taxonomic tree, not a sequence phylogeny).

## Outputs

- `itol_output_reference/` — the Newick tree + iTOL annotation files (substrate-composition
  tracks).
- `protein_3group_violin.png`, `protein_cellularity_stats.csv` — the cellularity-group
  comparison figure and its statistics.

## Why it matters

The classifier is only useful if its predictions recapitulate known biology at scale.
These trends (e.g. enrichment of protein-substrate MTases in more complex eukaryotes,
lineage-specific substrate biases) are the macro-scale sanity check and a source of
biological hypotheses from the tree-of-life sweep.
