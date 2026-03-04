# Phylogenetics

Phylogenetic tree construction and visualization for SAM MTase substrate preferences across organisms.

## Scripts

| Script | Description |
|--------|-------------|
| `reference_tree.py` | Builds iTOL-compatible phylogenetic tree with substrate annotations. Parses taxonomic lineages, counts SAMs per organism (≥20 minimum), and generates gradient/bar annotation files. |
| `mycobacterium_tree.py` | Mycobacterium-specific phylogenetic tree construction and annotation. |
| `grabbing_phylo_tree.py` | Retrieves phylogenetic tree data from external sources. |
| `grabbing_phylo_subtree.py` | Extracts and processes phylogenetic subtrees. |
| `grabbing_phylo_subtree_percent_SAM.py` | Subtree extraction with SAM percentage annotations. |
| `make_graphlan_annotations.py` | Generates GraPhlAn-format annotation files for circular tree visualization. |
| `uniprot_API.py` | UniProt API interface for retrieving protein and proteome metadata. |

## Dependencies

- [iTOL](https://itol.embl.de/) — Interactive Tree of Life for tree visualization
- [GraPhlAn](https://github.com/biobakery/graphlan) — Circular phylogenetic tree visualization (reference only; install separately)

## Output

- Newick format phylogenetic trees
- iTOL annotation files (DATASET_GRADIENT, DATASET_SIMPLEBAR)
- GraPhlAn XML and annotation files

## Color Scheme

| Substrate | Color |
|-----------|-------|
| RNA       | Red   |
| DNA       | Purple |
| Protein   | Green |
| Other     | Orange |
