#!/usr/bin/env python3
"""
SAM MTase Substrate Expression Analysis - Visualization Script

Generates publication-quality figures from the analysis results.
Run this after sam_mtase_substrate_expression_analysis.py has been executed.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
import sys

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# Plot styling
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 10

# Color palette for substrate types
SUBSTRATE_COLORS = {
    'RNA': '#E74C3C',      # Red
    'DNA': '#3498DB',      # Blue
    'Protein': '#2ECC71',  # Green
    'Other': '#9B59B6'     # Purple
}


def load_analysis_results(results_dir: Path) -> dict:
    """Load all analysis result files."""
    logger.info(f"Loading results from: {results_dir}")

    results = {}

    files_to_load = {
        'heatmap_mean': 'expression_heatmap_mean.csv',
        'heatmap_median': 'expression_heatmap_median.csv',
        'expression_summary': 'expression_summary_by_substrate_tissue.csv',
        'tau_summary': 'tau_summary_by_substrate.csv',
        'protein_tau': 'protein_tau_and_max_tissue.csv',
        'max_tissue_crosstab': 'max_expression_tissue_by_substrate.csv',
        'kruskal_by_tissue': 'kruskal_wallis_by_tissue.csv',
        'enrichment_results': 'enrichment_high_expressors.csv',
        'long_df': 'expression_long_format.csv',
        'coexpression_summary': 'coexpression_summary.csv',
    }

    for key, filename in files_to_load.items():
        filepath = results_dir / filename
        if filepath.exists():
            if key in ['heatmap_mean', 'heatmap_median', 'max_tissue_crosstab']:
                results[key] = pd.read_csv(filepath, index_col=0)
            elif key == 'tau_summary':
                results[key] = pd.read_csv(filepath, index_col=0)
            else:
                results[key] = pd.read_csv(filepath)
            logger.info(f"  Loaded: {filename}")
        else:
            logger.warning(f"  Missing: {filename}")

    return results


# =============================================================================
# FIGURE 1: Expression Heatmap
# =============================================================================

def plot_expression_heatmap(results: dict, output_dir: Path):
    """Create heatmap of mean expression by substrate type across tissues."""
    logger.info("Creating expression heatmap...")

    if 'heatmap_mean' not in results:
        logger.warning("Heatmap data not found, skipping...")
        return

    heatmap_data = results['heatmap_mean']

    # Log-transform for better visualization (add pseudocount)
    heatmap_log = np.log10(heatmap_data + 1)

    fig, ax = plt.subplots(figsize=(16, 4))

    sns.heatmap(
        heatmap_log,
        cmap='YlOrRd',
        ax=ax,
        cbar_kws={'label': 'log10(nTPM + 1)'},
        linewidths=0.5
    )

    ax.set_title('Mean Expression by Substrate Type Across Tissues', fontsize=14, fontweight='bold')
    ax.set_xlabel('Tissue', fontsize=11)
    ax.set_ylabel('Predicted Substrate Type', fontsize=11)

    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.yticks(fontsize=10)

    plt.tight_layout()
    plt.savefig(output_dir / 'fig1_expression_heatmap.png', bbox_inches='tight')
    plt.savefig(output_dir / 'fig1_expression_heatmap.pdf', bbox_inches='tight')
    plt.close()

    logger.info("  Saved: fig1_expression_heatmap.png/pdf")


# =============================================================================
# FIGURE 2: Tissue Specificity (Tau) Distributions
# =============================================================================

def plot_tau_distributions(results: dict, output_dir: Path):
    """Create violin/box plots of Tau distributions by substrate type."""
    logger.info("Creating Tau distribution plots...")

    if 'protein_tau' not in results:
        logger.warning("Tau data not found, skipping...")
        return

    df = results['protein_tau']

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Plot 1: Violin plot
    ax1 = axes[0]
    order = ['RNA', 'DNA', 'Protein', 'Other']
    available_order = [s for s in order if s in df['Predicted_Substrate'].unique()]

    palette = [SUBSTRATE_COLORS.get(s, '#888888') for s in available_order]

    sns.violinplot(
        data=df,
        x='Predicted_Substrate',
        y='Tau',
        order=available_order,
        palette=palette,
        ax=ax1,
        inner='box'
    )

    ax1.set_title('Tissue Specificity (Tau) by Substrate Type', fontsize=12, fontweight='bold')
    ax1.set_xlabel('Predicted Substrate Type', fontsize=11)
    ax1.set_ylabel('Tau Index\n(0 = ubiquitous, 1 = tissue-specific)', fontsize=10)
    ax1.set_ylim(0, 1)

    # Add sample sizes
    for i, substrate in enumerate(available_order):
        n = len(df[df['Predicted_Substrate'] == substrate])
        ax1.text(i, -0.08, f'n={n}', ha='center', fontsize=9, transform=ax1.get_xaxis_transform())

    # Plot 2: Box plot with individual points
    ax2 = axes[1]

    sns.boxplot(
        data=df,
        x='Predicted_Substrate',
        y='Tau',
        order=available_order,
        palette=palette,
        ax=ax2,
        width=0.5
    )

    sns.stripplot(
        data=df,
        x='Predicted_Substrate',
        y='Tau',
        order=available_order,
        color='black',
        alpha=0.3,
        size=3,
        ax=ax2
    )

    ax2.set_title('Tau Distribution with Individual Points', fontsize=12, fontweight='bold')
    ax2.set_xlabel('Predicted Substrate Type', fontsize=11)
    ax2.set_ylabel('Tau Index', fontsize=10)
    ax2.set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig(output_dir / 'fig2_tau_distributions.png', bbox_inches='tight')
    plt.savefig(output_dir / 'fig2_tau_distributions.pdf', bbox_inches='tight')
    plt.close()

    logger.info("  Saved: fig2_tau_distributions.png/pdf")


# =============================================================================
# FIGURE 3: Tissue of Maximum Expression
# =============================================================================

def plot_max_expression_tissue(results: dict, output_dir: Path):
    """Stacked bar chart showing tissue of max expression by substrate type."""
    logger.info("Creating max expression tissue plot...")

    if 'max_tissue_crosstab' not in results:
        logger.warning("Max tissue data not found, skipping...")
        return

    crosstab = results['max_tissue_crosstab']

    # Sort tissues by total frequency
    tissue_totals = crosstab.sum(axis=0).sort_values(ascending=False)
    top_tissues = tissue_totals.head(15).index.tolist()  # Top 15 tissues

    # Subset and add "Other" category
    crosstab_plot = crosstab[top_tissues].copy()
    other_tissues = [c for c in crosstab.columns if c not in top_tissues]
    if other_tissues:
        crosstab_plot['Other tissues'] = crosstab[other_tissues].sum(axis=1)

    fig, ax = plt.subplots(figsize=(14, 6))

    crosstab_plot.T.plot(
        kind='bar',
        stacked=True,
        ax=ax,
        color=[SUBSTRATE_COLORS.get(s, '#888888') for s in crosstab_plot.index],
        edgecolor='white',
        linewidth=0.5
    )

    ax.set_title('Tissue of Maximum Expression by Substrate Type', fontsize=14, fontweight='bold')
    ax.set_xlabel('Tissue', fontsize=11)
    ax.set_ylabel('Percentage of Proteins', fontsize=11)
    ax.legend(title='Substrate Type', bbox_to_anchor=(1.02, 1), loc='upper left')

    plt.xticks(rotation=45, ha='right', fontsize=9)

    plt.tight_layout()
    plt.savefig(output_dir / 'fig3_max_expression_tissue.png', bbox_inches='tight')
    plt.savefig(output_dir / 'fig3_max_expression_tissue.pdf', bbox_inches='tight')
    plt.close()

    logger.info("  Saved: fig3_max_expression_tissue.png/pdf")


# =============================================================================
# FIGURE 4: Expression Box Plots for Selected Tissues
# =============================================================================

def plot_expression_boxplots(results: dict, output_dir: Path):
    """Box plots of expression by substrate type for selected tissues."""
    logger.info("Creating expression box plots...")

    if 'long_df' not in results:
        logger.warning("Long format data not found, skipping...")
        return

    df = results['long_df']

    # Select interesting tissues (brain, liver, testis, etc.)
    interesting_tissues = [
        'cerebral cortex', 'liver', 'testis', 'heart muscle',
        'skeletal muscle', 'bone marrow', 'thymus', 'kidney'
    ]

    available_tissues = [t for t in interesting_tissues if t in df['Tissue'].unique()]

    if len(available_tissues) == 0:
        # Fallback: use most variable tissues
        tissue_var = df.groupby('Tissue')['nTPM'].var().sort_values(ascending=False)
        available_tissues = tissue_var.head(8).index.tolist()

    n_tissues = len(available_tissues)
    n_cols = 4
    n_rows = (n_tissues + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4 * n_rows))
    axes = axes.flatten() if n_rows > 1 or n_cols > 1 else [axes]

    order = ['RNA', 'DNA', 'Protein', 'Other']
    available_order = [s for s in order if s in df['Predicted_Substrate'].unique()]
    palette = [SUBSTRATE_COLORS.get(s, '#888888') for s in available_order]

    for i, tissue in enumerate(available_tissues):
        ax = axes[i]
        tissue_data = df[df['Tissue'] == tissue]

        sns.boxplot(
            data=tissue_data,
            x='Predicted_Substrate',
            y='nTPM',
            order=available_order,
            palette=palette,
            ax=ax
        )

        ax.set_title(tissue.title(), fontsize=11, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('nTPM' if i % n_cols == 0 else '')
        ax.set_yscale('symlog', linthresh=1)  # Symlog for better visualization

    # Hide unused axes
    for i in range(len(available_tissues), len(axes)):
        axes[i].set_visible(False)

    fig.suptitle('Expression Distribution by Substrate Type Across Tissues',
                 fontsize=14, fontweight='bold', y=1.02)

    plt.tight_layout()
    plt.savefig(output_dir / 'fig4_expression_boxplots.png', bbox_inches='tight')
    plt.savefig(output_dir / 'fig4_expression_boxplots.pdf', bbox_inches='tight')
    plt.close()

    logger.info("  Saved: fig4_expression_boxplots.png/pdf")


# =============================================================================
# FIGURE 5: Statistical Significance Overview
# =============================================================================

def plot_statistical_significance(results: dict, output_dir: Path):
    """Heatmap showing statistical significance of enrichment tests."""
    logger.info("Creating statistical significance heatmap...")

    if 'enrichment_results' not in results:
        logger.warning("Enrichment results not found, skipping...")
        return

    enrich_df = results['enrichment_results']

    # Create pivot table of -log10(p-value) * sign(log(OR))
    enrich_df['signed_log_p'] = -np.log10(enrich_df['p_value'] + 1e-300) * np.sign(np.log(enrich_df['Odds_Ratio'] + 0.01))

    pivot = enrich_df.pivot(
        index='Substrate',
        columns='Tissue',
        values='signed_log_p'
    )

    # Order substrates
    order = ['RNA', 'DNA', 'Protein', 'Other']
    available_order = [s for s in order if s in pivot.index]
    pivot = pivot.reindex(available_order)

    fig, ax = plt.subplots(figsize=(18, 4))

    # Custom colormap: blue (depleted) - white - red (enriched)
    cmap = sns.diverging_palette(240, 10, as_cmap=True)

    # Clip extreme values for visualization
    vmax = min(10, pivot.abs().max().max())

    sns.heatmap(
        pivot,
        cmap=cmap,
        center=0,
        vmin=-vmax,
        vmax=vmax,
        ax=ax,
        cbar_kws={'label': '-log10(p) × sign(log(OR))\n(+: enriched, -: depleted)'},
        linewidths=0.5
    )

    ax.set_title('Enrichment of Substrate Types Among High Expressors (Top 25%)',
                 fontsize=14, fontweight='bold')
    ax.set_xlabel('Tissue', fontsize=11)
    ax.set_ylabel('Substrate Type', fontsize=11)

    plt.xticks(rotation=45, ha='right', fontsize=8)

    plt.tight_layout()
    plt.savefig(output_dir / 'fig5_enrichment_heatmap.png', bbox_inches='tight')
    plt.savefig(output_dir / 'fig5_enrichment_heatmap.pdf', bbox_inches='tight')
    plt.close()

    logger.info("  Saved: fig5_enrichment_heatmap.png/pdf")


# =============================================================================
# FIGURE 6: Co-expression Summary
# =============================================================================

def plot_coexpression_summary(results: dict, output_dir: Path):
    """Bar plot comparing within vs between group correlations."""
    logger.info("Creating co-expression summary plot...")

    if 'coexpression_summary' not in results:
        logger.warning("Co-expression summary not found, skipping...")
        return

    coexp_df = results['coexpression_summary']

    # Extract values
    within_vals = {}
    for _, row in coexp_df.iterrows():
        metric = row['Metric']
        value = row['Value']

        if 'within-group mean' in metric and metric.startswith('Overall'):
            within_vals['All (within)'] = value
        elif 'between-group mean' in metric:
            within_vals['All (between)'] = value
        elif 'within-group mean correlation' in metric:
            substrate = metric.split(' within')[0]
            within_vals[substrate] = value

    fig, ax = plt.subplots(figsize=(10, 6))

    # Prepare data for plotting
    categories = []
    values = []
    colors = []

    # Add substrate-specific within-group correlations
    for substrate in ['RNA', 'DNA', 'Protein', 'Other']:
        if substrate in within_vals:
            categories.append(f'{substrate}\n(within)')
            values.append(within_vals[substrate])
            colors.append(SUBSTRATE_COLORS.get(substrate, '#888888'))

    # Add overall comparison
    if 'All (within)' in within_vals:
        categories.append('All\n(within)')
        values.append(within_vals['All (within)'])
        colors.append('#333333')

    if 'All (between)' in within_vals:
        categories.append('All\n(between)')
        values.append(within_vals['All (between)'])
        colors.append('#999999')

    bars = ax.bar(categories, values, color=colors, edgecolor='black', linewidth=0.5)

    ax.set_title('Mean Expression Profile Correlation\nWithin vs Between Substrate Types',
                 fontsize=14, fontweight='bold')
    ax.set_ylabel('Mean Pearson Correlation', fontsize=11)
    ax.set_xlabel('Group', fontsize=11)

    # Add value labels
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f'{val:.3f}', ha='center', va='bottom', fontsize=9)

    ax.set_ylim(0, max(values) * 1.15)

    plt.tight_layout()
    plt.savefig(output_dir / 'fig6_coexpression_summary.png', bbox_inches='tight')
    plt.savefig(output_dir / 'fig6_coexpression_summary.pdf', bbox_inches='tight')
    plt.close()

    logger.info("  Saved: fig6_coexpression_summary.png/pdf")


# =============================================================================
# FIGURE 7: Summary Overview Panel
# =============================================================================

def plot_summary_panel(results: dict, output_dir: Path):
    """Create a multi-panel summary figure."""
    logger.info("Creating summary panel figure...")

    fig = plt.figure(figsize=(16, 12))

    # Panel A: Substrate distribution (pie chart)
    ax1 = fig.add_subplot(2, 3, 1)
    if 'protein_tau' in results:
        df = results['protein_tau']
        counts = df['Predicted_Substrate'].value_counts()
        order = ['RNA', 'DNA', 'Protein', 'Other']
        counts = counts.reindex([s for s in order if s in counts.index])
        colors = [SUBSTRATE_COLORS.get(s, '#888888') for s in counts.index]

        ax1.pie(counts.values, labels=counts.index, colors=colors,
                autopct='%1.1f%%', startangle=90)
        ax1.set_title('A. Substrate Type Distribution', fontsize=12, fontweight='bold')

    # Panel B: Tau by substrate (box plot)
    ax2 = fig.add_subplot(2, 3, 2)
    if 'protein_tau' in results:
        df = results['protein_tau']
        order = ['RNA', 'DNA', 'Protein', 'Other']
        available_order = [s for s in order if s in df['Predicted_Substrate'].unique()]
        palette = [SUBSTRATE_COLORS.get(s, '#888888') for s in available_order]

        sns.boxplot(data=df, x='Predicted_Substrate', y='Tau',
                    order=available_order, palette=palette, ax=ax2)
        ax2.set_title('B. Tissue Specificity (Tau)', fontsize=12, fontweight='bold')
        ax2.set_xlabel('')
        ax2.set_ylabel('Tau Index')
        ax2.set_ylim(0, 1)

    # Panel C: Mean expression heatmap (subset)
    ax3 = fig.add_subplot(2, 3, 3)
    if 'heatmap_mean' in results:
        heatmap = results['heatmap_mean']
        # Select subset of tissues
        interesting = ['cerebral cortex', 'liver', 'testis', 'heart muscle',
                       'skeletal muscle', 'kidney', 'thymus', 'bone marrow']
        available = [t for t in interesting if t in heatmap.columns]
        if len(available) < 4:
            available = heatmap.columns[:8].tolist()

        heatmap_sub = heatmap[available]
        heatmap_log = np.log10(heatmap_sub + 1)

        sns.heatmap(heatmap_log, cmap='YlOrRd', ax=ax3,
                    cbar_kws={'label': 'log10(nTPM+1)'})
        ax3.set_title('C. Mean Expression (Selected Tissues)', fontsize=12, fontweight='bold')
        ax3.set_xlabel('')
        plt.sca(ax3)
        plt.xticks(rotation=45, ha='right', fontsize=8)

    # Panel D: Significant tissues bar
    ax4 = fig.add_subplot(2, 3, 4)
    if 'kruskal_by_tissue' in results:
        kruskal = results['kruskal_by_tissue']
        kruskal_sig = kruskal[kruskal['significant']].sort_values('H_statistic', ascending=False).head(10)

        if len(kruskal_sig) > 0:
            ax4.barh(kruskal_sig['Tissue'], kruskal_sig['H_statistic'], color='steelblue')
            ax4.set_xlabel('Kruskal-Wallis H Statistic')
            ax4.set_title('D. Tissues with Significant\nExpression Differences', fontsize=12, fontweight='bold')
            ax4.invert_yaxis()
        else:
            ax4.text(0.5, 0.5, 'No significant\ntissues found', ha='center', va='center', transform=ax4.transAxes)
            ax4.set_title('D. Tissues with Significant Differences', fontsize=12, fontweight='bold')

    # Panel E: Top enrichments
    ax5 = fig.add_subplot(2, 3, 5)
    if 'enrichment_results' in results:
        enrich = results['enrichment_results']
        enriched = enrich[enrich['enriched']].nlargest(10, 'Odds_Ratio')

        if len(enriched) > 0:
            labels = enriched['Substrate'] + ' in ' + enriched['Tissue']
            colors = [SUBSTRATE_COLORS.get(s, '#888888') for s in enriched['Substrate']]
            ax5.barh(labels, enriched['Odds_Ratio'], color=colors)
            ax5.set_xlabel('Odds Ratio')
            ax5.set_title('E. Top Substrate Enrichments\n(High Expressors)', fontsize=12, fontweight='bold')
            ax5.invert_yaxis()
            ax5.axvline(x=1, color='black', linestyle='--', linewidth=0.5)
        else:
            ax5.text(0.5, 0.5, 'No significant\nenrichments found', ha='center', va='center', transform=ax5.transAxes)
            ax5.set_title('E. Top Substrate Enrichments', fontsize=12, fontweight='bold')

    # Panel F: Co-expression comparison
    ax6 = fig.add_subplot(2, 3, 6)
    if 'coexpression_summary' in results:
        coexp = results['coexpression_summary']
        within = coexp[coexp['Metric'] == 'Overall within-group mean correlation']['Value'].values
        between = coexp[coexp['Metric'] == 'Overall between-group mean correlation']['Value'].values

        if len(within) > 0 and len(between) > 0:
            bars = ax6.bar(['Within-group', 'Between-group'], [within[0], between[0]],
                           color=['#2ECC71', '#E74C3C'], edgecolor='black')
            ax6.set_ylabel('Mean Correlation')
            ax6.set_title('F. Expression Profile Correlation', fontsize=12, fontweight='bold')

            for bar, val in zip(bars, [within[0], between[0]]):
                ax6.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                         f'{val:.3f}', ha='center', va='bottom', fontsize=10)

    plt.tight_layout()
    plt.savefig(output_dir / 'fig7_summary_panel.png', bbox_inches='tight')
    plt.savefig(output_dir / 'fig7_summary_panel.pdf', bbox_inches='tight')
    plt.close()

    logger.info("  Saved: fig7_summary_panel.png/pdf")


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    """Main execution function."""
    logger.info("=" * 60)
    logger.info("SAM MTase Substrate Expression - Visualization")
    logger.info("=" * 60)

    # Define paths
    script_dir = Path(__file__).parent
    results_dir = script_dir / "substrate_expression_analysis"
    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    if not results_dir.exists():
        logger.error(f"Results directory not found: {results_dir}")
        logger.error("Please run sam_mtase_substrate_expression_analysis.py first!")
        return

    # Load results
    results = load_analysis_results(results_dir)

    if not results:
        logger.error("No results loaded!")
        return

    # Generate all figures
    logger.info("\nGenerating figures...")

    plot_expression_heatmap(results, figures_dir)
    plot_tau_distributions(results, figures_dir)
    plot_max_expression_tissue(results, figures_dir)
    plot_expression_boxplots(results, figures_dir)
    plot_statistical_significance(results, figures_dir)
    plot_coexpression_summary(results, figures_dir)
    plot_summary_panel(results, figures_dir)

    logger.info("=" * 60)
    logger.info(f"All figures saved to: {figures_dir}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
