# SAM MTase Substrate Expression Analysis

This directory contains analysis results examining the relationship between predicted SAM methyltransferase (MTase) substrate types and tissue expression patterns in humans.

## Overview

The analysis was performed by two scripts:
1. **sam_mtase_substrate_expression_analysis.py** - Performs all statistical analyses and generates CSV data files
2. **sam_mtase_substrate_expression_plots.py** - Generates publication-quality figures from the analysis results

### Input Data
- **Expression data**: `sam_mtase_tissue_expression_matrix.csv` - Tissue expression matrix (nTPM values) for human SAM MTases across 51 tissues
- **Predictions**: Human substrate type predictions from the substrate type predictor (probabilities for RNA, DNA, Protein, Other)

### Dataset Summary
- **374 proteins** analyzed (proteins with both expression data and substrate predictions)
- **51 tissues** examined
- **4 predicted substrate types**: RNA (27.8%), DNA (5.1%), Protein (39.6%), Other (27.5%)

---

## CSV Files

### 1. `expression_heatmap_mean.csv`
**Description**: Matrix of mean expression values (nTPM) for each substrate type across all tissues.

**Structure**:
- Rows: Predicted substrate types (DNA, Other, Protein, RNA)
- Columns: Tissue names (51 tissues)
- Values: Mean nTPM expression

**How it was calculated**: For each substrate type and tissue combination, the arithmetic mean of nTPM values was computed across all proteins of that substrate type.

---

### 2. `expression_heatmap_median.csv`
**Description**: Matrix of median expression values (nTPM) for each substrate type across all tissues.

**Structure**: Same as `expression_heatmap_mean.csv`

**How it was calculated**: For each substrate type and tissue combination, the median nTPM value was computed across all proteins of that substrate type.

---

### 3. `expression_summary_by_substrate_tissue.csv`
**Description**: Comprehensive summary statistics for expression by substrate type and tissue.

**Columns**:
- `Predicted_Substrate`: The predicted substrate type
- `Tissue`: Tissue name
- `mean`: Mean nTPM expression
- `median`: Median nTPM expression
- `std`: Standard deviation
- `count`: Number of proteins
- `sem`: Standard error of the mean (std / sqrt(n))

**How it was calculated**: Expression data was melted to long format, then grouped by substrate type and tissue to compute summary statistics.

---

### 4. `tau_summary_by_substrate.csv`
**Description**: Summary of tissue specificity (Tau index) by substrate type.

**Columns**:
- Index: Substrate type
- `mean`: Mean Tau value
- `median`: Median Tau value
- `std`: Standard deviation
- `count`: Number of proteins

**How it was calculated**:
The Tau index measures tissue specificity on a scale from 0 (ubiquitous expression) to 1 (tissue-specific expression):

```
Tau = sum(1 - x_hat_i) / (n - 1)
where x_hat_i = x_i / max(x)
```

For each protein, Tau was calculated from its expression profile across all tissues. Then Tau values were grouped by substrate type to compute summary statistics.

---

### 5. `protein_tau_and_max_tissue.csv`
**Description**: Per-protein information including Tau values and tissue of maximum expression.

**Columns**:
- `Entry`: UniProt accession
- `Protein names`: Protein name(s)
- `Predicted_Substrate`: Predicted substrate type (highest probability class)
- `Prediction_Confidence`: Probability of the predicted class
- `Tau`: Tissue specificity index (0-1)
- `Max_Expression_Tissue`: Tissue with highest expression for this protein
- `Max_Expression_Value`: nTPM value in the max-expressing tissue

**How it was calculated**: For each protein, Tau was computed from its expression vector. The tissue with the maximum expression value was identified along with that value.

---

### 6. `max_expression_tissue_by_substrate.csv`
**Description**: Cross-tabulation showing the percentage of proteins in each substrate type that have their maximum expression in each tissue.

**Structure**:
- Rows: Predicted substrate types
- Columns: Tissue names
- Values: Percentage of proteins (within each substrate type) that have max expression in that tissue

**How it was calculated**: Created using `pd.crosstab()` with row-wise normalization (`normalize='index'`), then multiplied by 100 to convert to percentages.

---

### 7. `kruskal_wallis_by_tissue.csv`
**Description**: Results of Kruskal-Wallis tests examining whether expression levels differ significantly between substrate types within each tissue.

**Columns**:
- `Tissue`: Tissue name
- `H_statistic`: Kruskal-Wallis H statistic
- `p_value`: Raw p-value
- `p_adjusted`: Bonferroni-corrected p-value (multiplied by number of tests)
- `significant`: Boolean indicating if p_adjusted < 0.05

**How it was calculated**: For each tissue, a Kruskal-Wallis H-test (non-parametric ANOVA) was performed comparing expression values across all substrate types. Multiple testing correction was applied using Bonferroni correction.

---

### 8. `enrichment_high_expressors.csv`
**Description**: Fisher's exact test results examining whether each substrate type is enriched or depleted among high expressors in each tissue.

**Columns**:
- `Tissue`: Tissue name
- `Substrate`: Substrate type being tested
- `High_Substrate`: Count of high expressors of this substrate type
- `High_Other`: Count of high expressors of other substrate types
- `Low_Substrate`: Count of low expressors of this substrate type
- `Low_Other`: Count of low expressors of other substrate types
- `Odds_Ratio`: Fisher's exact test odds ratio
- `p_value`: Raw p-value
- `p_adjusted`: Bonferroni-corrected p-value
- `significant`: Boolean (p_adjusted < 0.05)
- `enriched`: Boolean (Odds_Ratio > 1 AND significant)
- `depleted`: Boolean (Odds_Ratio < 1 AND significant)

**How it was calculated**:
1. For each tissue, proteins were classified as "high expressors" (top 25%, >= 75th percentile) or not
2. A 2x2 contingency table was constructed for each substrate type (substrate vs. other, high vs. low)
3. Fisher's exact test was performed to test for enrichment/depletion
4. Bonferroni correction was applied for multiple testing

---

### 9. `pairwise_mannwhitney_tests.csv`
**Description**: Pairwise Mann-Whitney U test results comparing expression between each pair of substrate types within each tissue.

**Columns**:
- `Tissue`: Tissue name
- `Comparison`: Comparison label (e.g., "RNA vs DNA")
- `Substrate_1`, `Substrate_2`: The two substrate types being compared
- `Mean_1`, `Mean_2`: Mean expression values for each group
- `U_statistic`: Mann-Whitney U statistic
- `p_value`: Raw p-value (two-sided test)
- `p_adjusted`: Bonferroni-corrected p-value
- `significant`: Boolean (p_adjusted < 0.05)

**How it was calculated**: For each tissue, all pairwise combinations of substrate types were compared using two-sided Mann-Whitney U tests. Bonferroni correction was applied across all tests.

---

### 10. `coexpression_summary.csv`
**Description**: Summary of co-expression analysis comparing expression profile correlations within vs. between substrate type groups.

**Columns**:
- `Metric`: Description of the metric
- `Value`: Numeric value

**Metrics included**:
- `Overall within-group mean correlation`: Mean Pearson correlation between proteins of the SAME substrate type
- `Overall between-group mean correlation`: Mean Pearson correlation between proteins of DIFFERENT substrate types
- `Within-group pairs count`: Number of within-group protein pairs analyzed
- `Between-group pairs count`: Number of between-group protein pairs analyzed
- `{Substrate} within-group mean correlation`: Mean correlation within each specific substrate type
- `{Substrate} within-group pair count`: Number of pairs within each substrate type

**How it was calculated**:
1. A pairwise Pearson correlation matrix was computed for all protein expression profiles
2. Correlations were separated into "within-group" (same substrate type) and "between-group" (different substrate types)
3. Summary statistics were computed for each category
4. A Mann-Whitney U test (one-sided, alternative='greater') was performed to test if within-group correlations are significantly higher

---

### 11. `expression_long_format.csv`
**Description**: Full expression data in long/tidy format for custom analyses or plotting.

**Columns**:
- `Entry`: UniProt accession
- `Protein names`: Protein name(s)
- `Predicted_Substrate`: Predicted substrate type
- `Prediction_Confidence`: Prediction probability
- `Tissue`: Tissue name
- `nTPM`: Expression value (normalized transcripts per million)

**How it was calculated**: The merged expression matrix was melted from wide format (proteins × tissues) to long format (one row per protein-tissue combination).

---

### 12. `analysis_summary_report.txt`
**Description**: Human-readable text summary of key findings from the analysis.

---

## Figures

All figures are provided in both PNG (300 dpi) and PDF formats in the `figures/` subdirectory.

### Figure 1: `fig1_expression_heatmap.png/pdf`
**Description**: Heatmap showing mean expression by substrate type across all tissues.

**Visualization details**:
- Y-axis: Predicted substrate types (DNA, Other, Protein, RNA)
- X-axis: Tissues (51)
- Color scale: log10(nTPM + 1), YlOrRd colormap
- The log transformation improves visualization of the wide dynamic range of expression values

**Data source**: `expression_heatmap_mean.csv`

---

### Figure 2: `fig2_tau_distributions.png/pdf`
**Description**: Two-panel figure showing tissue specificity (Tau) distributions by substrate type.

**Panel 1 (Left)**: Violin plot with embedded box plot showing Tau distribution for each substrate type
**Panel 2 (Right)**: Box plot with overlaid individual data points (strip plot)

**Interpretation**: Higher Tau values indicate more tissue-specific expression. The figure reveals whether certain substrate types tend to have more ubiquitous or tissue-specific expression patterns.

**Data source**: `protein_tau_and_max_tissue.csv`

---

### Figure 3: `fig3_max_expression_tissue.png/pdf`
**Description**: Stacked bar chart showing which tissues have the highest expression for proteins of each substrate type.

**Visualization details**:
- X-axis: Top 15 most common tissues of max expression (plus "Other tissues")
- Y-axis: Percentage of proteins
- Colors: Stacked by substrate type (RNA=red, DNA=blue, Protein=green, Other=purple)

**Data source**: `max_expression_tissue_by_substrate.csv`

---

### Figure 4: `fig4_expression_boxplots.png/pdf`
**Description**: Multi-panel figure showing expression distributions by substrate type for selected tissues of biological interest.

**Tissues shown** (if available): cerebral cortex, liver, testis, heart muscle, skeletal muscle, bone marrow, thymus, kidney

**Visualization details**:
- Each panel shows one tissue
- X-axis: Substrate type
- Y-axis: nTPM (symlog scale with linear threshold at 1)
- Colors match substrate type color scheme

**Data source**: `expression_long_format.csv`

---

### Figure 5: `fig5_enrichment_heatmap.png/pdf`
**Description**: Heatmap showing statistical significance of substrate type enrichment among high expressors.

**Visualization details**:
- Y-axis: Substrate types
- X-axis: Tissues
- Color scale: Diverging colormap representing -log10(p) × sign(log(OR))
  - Red: Substrate type enriched among top 25% expressors
  - Blue: Substrate type depleted among top 25% expressors
  - White: No significant enrichment/depletion
- Color intensity reflects statistical significance

**Data source**: `enrichment_high_expressors.csv`

---

### Figure 6: `fig6_coexpression_summary.png/pdf`
**Description**: Bar chart comparing mean expression profile correlations within vs. between substrate type groups.

**Bars shown**:
- Individual substrate types (within-group correlation)
- "All (within)": Overall mean within-group correlation
- "All (between)": Overall mean between-group correlation

**Interpretation**: If proteins of the same substrate type have similar biological functions, we expect within-group correlations to be higher than between-group correlations.

**Data source**: `coexpression_summary.csv`

---

### Figure 7: `fig7_summary_panel.png/pdf`
**Description**: Six-panel summary figure providing an overview of all key findings.

**Panel A**: Pie chart of substrate type distribution
**Panel B**: Box plot of Tau by substrate type
**Panel C**: Expression heatmap for selected tissues
**Panel D**: Horizontal bar chart of tissues with significant expression differences (by Kruskal-Wallis H statistic)
**Panel E**: Horizontal bar chart of top substrate enrichments among high expressors (by odds ratio)
**Panel F**: Bar chart comparing within-group vs. between-group expression correlations

**Data sources**: Multiple CSV files

---

## Statistical Methods Summary

| Analysis | Test | Correction |
|----------|------|------------|
| Expression differences by substrate type per tissue | Kruskal-Wallis H | Bonferroni |
| Enrichment among high expressors | Fisher's exact (two-sided) | Bonferroni |
| Pairwise substrate comparisons | Mann-Whitney U (two-sided) | Bonferroni |
| Within vs. between-group correlations | Mann-Whitney U (one-sided) | None |
| Tau differences across substrate types | Kruskal-Wallis H | None (single test) |

---

## Color Scheme

Consistent colors are used throughout for substrate types:
- **RNA**: Red (#E74C3C)
- **DNA**: Blue (#3498DB)
- **Protein**: Green (#2ECC71)
- **Other**: Purple (#9B59B6)

---

## Reproducibility

To regenerate these results:
```bash
# From the graphlan directory:
python sam_mtase_substrate_expression_analysis.py  # Generate CSV files
python sam_mtase_substrate_expression_plots.py     # Generate figures
```

**Requirements**: pandas, numpy, scipy, matplotlib, seaborn
