#!/usr/bin/env python3
"""
Generate two heatmaps and a CSV for SAM MTase substrate expression across tissues:

1. Total expression: sum of nTPM for all enzymes in each substrate type per tissue
2. Percent makeup: each substrate type's share of total expression per tissue
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Plot styling
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

SUBSTRATE_ORDER = ['RNA', 'DNA', 'Protein', 'Other']

script_dir = Path(__file__).parent
results_dir = script_dir / "substrate_expression_analysis"
figures_dir = results_dir / "figures"
figures_dir.mkdir(parents=True, exist_ok=True)

# Load long-format expression data
long_df = pd.read_csv(results_dir / "expression_long_format.csv")

# --- Total summed expression per substrate type per tissue ---
total_expr = (
    long_df.groupby(['Predicted_Substrate', 'Tissue'])['nTPM']
    .sum()
    .reset_index()
    .rename(columns={'nTPM': 'Total_nTPM'})
)

# Pivot: substrate types as rows, tissues as columns
total_pivot = total_expr.pivot(
    index='Predicted_Substrate', columns='Tissue', values='Total_nTPM'
).reindex([s for s in SUBSTRATE_ORDER if s in total_expr['Predicted_Substrate'].unique()])

# Save CSV
total_pivot.to_csv(results_dir / "total_expression_by_substrate_tissue.csv")
print(f"Saved: {results_dir / 'total_expression_by_substrate_tissue.csv'}")

# --- Percent makeup: normalize each tissue column to sum to 100% ---
pct_pivot = total_pivot.div(total_pivot.sum(axis=0), axis=1) * 100

pct_pivot.to_csv(results_dir / "percent_expression_by_substrate_tissue.csv")
print(f"Saved: {results_dir / 'percent_expression_by_substrate_tissue.csv'}")

# --- Figure 1: Total expression heatmap ---
fig, ax = plt.subplots(figsize=(18, 4))

sns.heatmap(
    np.log10(total_pivot + 1),
    cmap='YlOrRd',
    ax=ax,
    cbar_kws={'label': 'log10(Total nTPM + 1)'},
    linewidths=0.5,
    annot=total_pivot.round(0).astype(int),
    fmt='d',
    annot_kws={'fontsize': 7},
)

ax.set_title('Total Expression by Substrate Type Across Tissues', fontsize=14, fontweight='bold')
ax.set_xlabel('Tissue', fontsize=11)
ax.set_ylabel('Substrate Type', fontsize=11)
plt.xticks(rotation=45, ha='right', fontsize=8)
plt.yticks(fontsize=10)
plt.tight_layout()

plt.savefig(figures_dir / 'total_expression_heatmap.png', bbox_inches='tight')
plt.savefig(figures_dir / 'total_expression_heatmap.pdf', bbox_inches='tight')
plt.close()
print(f"Saved: {figures_dir / 'total_expression_heatmap.png/pdf'}")

# --- Figure 2: Percent makeup heatmap ---
fig, ax = plt.subplots(figsize=(18, 4))

sns.heatmap(
    pct_pivot,
    cmap='coolwarm',
    center=25,  # 4 substrate types -> 25% would be even
    ax=ax,
    cbar_kws={'label': '% of Total Tissue Expression'},
    linewidths=0.5,
    annot=pct_pivot.round(1),
    fmt='.1f',
    annot_kws={'fontsize': 7},
    vmin=0,
    vmax=50,
)

ax.set_title('Percent Makeup of Expression by Substrate Type Across Tissues',
             fontsize=14, fontweight='bold')
ax.set_xlabel('Tissue', fontsize=11)
ax.set_ylabel('Substrate Type', fontsize=11)
plt.xticks(rotation=45, ha='right', fontsize=8)
plt.yticks(fontsize=10)
plt.tight_layout()

plt.savefig(figures_dir / 'percent_expression_heatmap.png', bbox_inches='tight')
plt.savefig(figures_dir / 'percent_expression_heatmap.pdf', bbox_inches='tight')
plt.close()
print(f"Saved: {figures_dir / 'percent_expression_heatmap.png/pdf'}")

# Color palette for substrate types
SUBSTRATE_COLORS = {
    'RNA': '#E74C3C',
    'DNA': '#3498DB',
    'Protein': '#2ECC71',
    'Other': '#9B59B6',
}

# --- Figure 3: Total expression stacked bar (sorted by Protein total) ---
total_sorted_by_protein = total_pivot.T.sort_values('Protein', ascending=False)

fig, ax = plt.subplots(figsize=(18, 6))
bottom = np.zeros(len(total_sorted_by_protein))

for substrate in SUBSTRATE_ORDER:
    if substrate in total_sorted_by_protein.columns:
        vals = total_sorted_by_protein[substrate].values
        ax.bar(total_sorted_by_protein.index, vals, bottom=bottom,
               label=substrate, color=SUBSTRATE_COLORS[substrate], edgecolor='white', linewidth=0.3)
        bottom += vals

ax.set_title('Total Expression by Substrate Type (Sorted by Protein Expression)',
             fontsize=14, fontweight='bold')
ax.set_xlabel('Tissue', fontsize=11)
ax.set_ylabel('Total nTPM', fontsize=11)
ax.legend(title='Substrate Type')
plt.xticks(rotation=45, ha='right', fontsize=8)
plt.tight_layout()

plt.savefig(figures_dir / 'total_expression_barplot.png', bbox_inches='tight')
plt.savefig(figures_dir / 'total_expression_barplot.pdf', bbox_inches='tight')
plt.close()
print(f"Saved: {figures_dir / 'total_expression_barplot.png/pdf'}")

# --- Figure 4: Percent makeup stacked bar (sorted by %Protein) ---
pct_sorted_by_protein = pct_pivot.T.sort_values('Protein', ascending=False)

fig, ax = plt.subplots(figsize=(18, 6))
bottom = np.zeros(len(pct_sorted_by_protein))

for substrate in SUBSTRATE_ORDER:
    if substrate in pct_sorted_by_protein.columns:
        vals = pct_sorted_by_protein[substrate].values
        ax.bar(pct_sorted_by_protein.index, vals, bottom=bottom,
               label=substrate, color=SUBSTRATE_COLORS[substrate], edgecolor='white', linewidth=0.3)
        bottom += vals

ax.set_title('Percent Makeup of Expression by Substrate Type (Sorted by % Protein)',
             fontsize=14, fontweight='bold')
ax.set_xlabel('Tissue', fontsize=11)
ax.set_ylabel('% of Total Expression', fontsize=11)
ax.set_ylim(0, 100)
ax.legend(title='Substrate Type')
plt.xticks(rotation=45, ha='right', fontsize=8)
plt.tight_layout()

plt.savefig(figures_dir / 'percent_expression_barplot.png', bbox_inches='tight')
plt.savefig(figures_dir / 'percent_expression_barplot.pdf', bbox_inches='tight')
plt.close()
print(f"Saved: {figures_dir / 'percent_expression_barplot.png/pdf'}")

print("\nDone!")
