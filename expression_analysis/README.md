# Expression Analysis

Human tissue expression analysis for SAM-dependent methyltransferases, examining relationships between predicted substrate type and expression patterns.

## Scripts

| Script | Description |
|--------|-------------|
| `sam_mtase_substrate_expression_analysis.py` | Comprehensive analysis: expression profiles by substrate type, tissue specificity (Tau index), statistical enrichment (Kruskal-Wallis, Fisher's exact), co-expression analysis, and pairwise Mann-Whitney U comparisons. |
| `sam_mtase_substrate_expression_plots.py` | Generates publication-quality expression figures (heatmaps, boxplots, bar charts). |
| `sam_mtase_tissue_expression.py` | Tissue-level expression processing and matrix generation. |
| `substrate_expression_totals.py` | Computes aggregate expression totals by substrate type across tissues. |

## Output

### Data Files
- Expression summary statistics by substrate type and tissue
- Tau index (tissue specificity: 0=ubiquitous, 1=tissue-specific) per protein
- Statistical test results (Kruskal-Wallis, Mann-Whitney, Fisher's exact)
- Co-expression correlation matrices

### Figures
Generated figures are saved to `../figures/expression/`:
- `fig1` — Expression heatmap
- `fig2` — Tau distributions by substrate type
- `fig3` — Max expression tissue assignments
- `fig4` — Expression boxplots
- `fig5` — Enrichment heatmap
- `fig6` — Co-expression summary
- `fig7` — Summary panel

## Usage

```bash
# Run full expression analysis
python sam_mtase_substrate_expression_analysis.py

# Generate figures
python sam_mtase_substrate_expression_plots.py
```
