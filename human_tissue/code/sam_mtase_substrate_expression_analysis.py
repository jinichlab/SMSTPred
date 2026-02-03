#!/usr/bin/env python3
"""
SAM MTase Substrate Type Expression Analysis

This script analyzes the relationship between predicted substrate types
(RNA, DNA, Protein, Other) and tissue expression patterns for human SAM MTases.

Analyses performed:
1. Expression Profile Comparison by Substrate Type
2. Tissue Specificity Analysis (Tau index)
3. Statistical Enrichment Tests
4. Co-expression Analysis
5. Biological Context Summary
"""

import pandas as pd
import numpy as np
import logging
import sys
from pathlib import Path
from scipy import stats
from typing import Tuple, Dict, List
import warnings
warnings.filterwarnings('ignore')

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)


# =============================================================================
# DATA LOADING
# =============================================================================

def load_expression_matrix(filepath: str) -> pd.DataFrame:
    """Load the expression matrix (proteins x tissues)."""
    logger.info(f"Loading expression matrix from: {filepath}")
    df = pd.read_csv(filepath, index_col='Entry')

    # Separate protein names from tissue columns
    tissue_cols = [c for c in df.columns if c != 'Protein names']

    logger.info(f"Loaded {len(df)} proteins x {len(tissue_cols)} tissues")
    return df


def load_predictions(filepath: str) -> pd.DataFrame:
    """Load substrate type predictions."""
    logger.info(f"Loading predictions from: {filepath}")
    df = pd.read_csv(filepath)

    # Columns of interest: Entry, Other, Protein, RNA, DNA
    pred_cols = ['Entry', 'Protein names', 'Other', 'Protein', 'RNA', 'DNA']
    available_cols = [c for c in pred_cols if c in df.columns]
    df = df[available_cols].copy()

    # Determine predicted substrate type (highest probability)
    substrate_cols = ['Other', 'Protein', 'RNA', 'DNA']
    df['Predicted_Substrate'] = df[substrate_cols].idxmax(axis=1)

    # Also store the confidence (max probability)
    df['Prediction_Confidence'] = df[substrate_cols].max(axis=1)

    logger.info(f"Loaded predictions for {len(df)} proteins")
    logger.info(f"Substrate distribution:\n{df['Predicted_Substrate'].value_counts()}")

    return df


def merge_expression_and_predictions(expression_df: pd.DataFrame,
                                      predictions_df: pd.DataFrame) -> pd.DataFrame:
    """Merge expression data with substrate predictions."""
    logger.info("Merging expression data with predictions...")

    # Reset index to merge on Entry
    expr_reset = expression_df.reset_index()

    # Merge
    merged = expr_reset.merge(
        predictions_df[['Entry', 'Predicted_Substrate', 'Prediction_Confidence',
                        'Other', 'Protein', 'RNA', 'DNA']],
        on='Entry',
        how='inner'
    )

    logger.info(f"Merged dataset: {len(merged)} proteins with both expression and predictions")

    return merged


# =============================================================================
# 1. EXPRESSION PROFILE COMPARISON BY SUBSTRATE TYPE
# =============================================================================

def calculate_expression_summary(merged_df: pd.DataFrame,
                                  tissue_cols: List[str]) -> pd.DataFrame:
    """Calculate mean/median expression per substrate type per tissue."""
    logger.info("Calculating expression summary by substrate type...")

    # Melt to long format
    id_vars = ['Entry', 'Protein names', 'Predicted_Substrate', 'Prediction_Confidence']
    long_df = merged_df.melt(
        id_vars=id_vars,
        value_vars=tissue_cols,
        var_name='Tissue',
        value_name='nTPM'
    )

    # Calculate summary statistics
    summary = long_df.groupby(['Predicted_Substrate', 'Tissue'])['nTPM'].agg([
        'mean', 'median', 'std', 'count',
        ('sem', lambda x: x.std() / np.sqrt(len(x)))
    ]).reset_index()

    return summary, long_df


def create_expression_heatmap_data(summary_df: pd.DataFrame) -> pd.DataFrame:
    """Create pivot table for heatmap visualization."""
    # Pivot: substrate types as rows, tissues as columns
    heatmap_mean = summary_df.pivot(
        index='Predicted_Substrate',
        columns='Tissue',
        values='mean'
    )

    heatmap_median = summary_df.pivot(
        index='Predicted_Substrate',
        columns='Tissue',
        values='median'
    )

    return heatmap_mean, heatmap_median


# =============================================================================
# 2. TISSUE SPECIFICITY ANALYSIS (TAU INDEX)
# =============================================================================

def calculate_tau(expression_vector: np.ndarray) -> float:
    """
    Calculate tissue specificity index (Tau).

    Tau = sum(1 - x_hat_i) / (n - 1)
    where x_hat_i = x_i / max(x)

    Tau near 0 = ubiquitous expression
    Tau near 1 = tissue-specific expression
    """
    x = np.array(expression_vector)

    # Handle all-zero case
    if x.max() == 0:
        return np.nan

    # Normalize by max
    x_hat = x / x.max()

    # Calculate Tau
    n = len(x)
    tau = np.sum(1 - x_hat) / (n - 1)

    return tau


def calculate_tau_for_all(merged_df: pd.DataFrame,
                          tissue_cols: List[str]) -> pd.DataFrame:
    """Calculate Tau for all proteins."""
    logger.info("Calculating tissue specificity (Tau) for all proteins...")

    tau_values = []
    max_tissue = []
    max_expression = []

    for idx, row in merged_df.iterrows():
        expr_vector = row[tissue_cols].values.astype(float)
        tau = calculate_tau(expr_vector)
        tau_values.append(tau)

        # Find tissue of maximum expression
        if np.nanmax(expr_vector) > 0:
            max_idx = np.nanargmax(expr_vector)
            max_tissue.append(tissue_cols[max_idx])
            max_expression.append(expr_vector[max_idx])
        else:
            max_tissue.append(None)
            max_expression.append(0)

    merged_df = merged_df.copy()
    merged_df['Tau'] = tau_values
    merged_df['Max_Expression_Tissue'] = max_tissue
    merged_df['Max_Expression_Value'] = max_expression

    logger.info(f"Tau calculated. Mean Tau = {np.nanmean(tau_values):.3f}")

    return merged_df


def compare_tau_by_substrate(merged_df: pd.DataFrame) -> pd.DataFrame:
    """Compare Tau distributions across substrate types."""
    logger.info("Comparing tissue specificity across substrate types...")

    tau_summary = merged_df.groupby('Predicted_Substrate')['Tau'].agg([
        'mean', 'median', 'std', 'count'
    ]).round(4)

    # Kruskal-Wallis test across all groups
    groups = [group['Tau'].dropna().values
              for name, group in merged_df.groupby('Predicted_Substrate')]

    if len(groups) >= 2 and all(len(g) > 0 for g in groups):
        stat, p_value = stats.kruskal(*groups)
        logger.info(f"Kruskal-Wallis test for Tau: H={stat:.3f}, p={p_value:.4e}")

    return tau_summary


def analyze_max_expression_tissue(merged_df: pd.DataFrame) -> pd.DataFrame:
    """Analyze which tissues show max expression for each substrate type."""
    logger.info("Analyzing tissue of maximum expression by substrate type...")

    # Cross-tabulation: substrate type vs max expression tissue
    crosstab = pd.crosstab(
        merged_df['Predicted_Substrate'],
        merged_df['Max_Expression_Tissue'],
        normalize='index'
    ) * 100  # Convert to percentage

    return crosstab


# =============================================================================
# 3. STATISTICAL ENRICHMENT TESTS
# =============================================================================

def test_expression_differences_per_tissue(long_df: pd.DataFrame) -> pd.DataFrame:
    """Test if expression differs by substrate type within each tissue."""
    logger.info("Testing expression differences per tissue (Kruskal-Wallis)...")

    results = []

    for tissue in long_df['Tissue'].unique():
        tissue_data = long_df[long_df['Tissue'] == tissue]

        groups = [group['nTPM'].values
                  for name, group in tissue_data.groupby('Predicted_Substrate')]

        if len(groups) >= 2 and all(len(g) > 0 for g in groups):
            stat, p_value = stats.kruskal(*groups)
            results.append({
                'Tissue': tissue,
                'H_statistic': stat,
                'p_value': p_value
            })

    results_df = pd.DataFrame(results)

    # Multiple testing correction (Bonferroni)
    results_df['p_adjusted'] = results_df['p_value'] * len(results_df)
    results_df['p_adjusted'] = results_df['p_adjusted'].clip(upper=1.0)
    results_df['significant'] = results_df['p_adjusted'] < 0.05

    results_df = results_df.sort_values('p_value')

    n_sig = results_df['significant'].sum()
    logger.info(f"Tissues with significant differences: {n_sig}/{len(results_df)}")

    return results_df


def test_enrichment_in_high_expressors(merged_df: pd.DataFrame,
                                        tissue_cols: List[str],
                                        percentile: int = 75) -> pd.DataFrame:
    """
    Test if substrate types are enriched among high expressors in each tissue.
    Uses Fisher's exact test.
    """
    logger.info(f"Testing substrate enrichment in top {100-percentile}% expressors per tissue...")

    results = []
    substrate_types = merged_df['Predicted_Substrate'].unique()

    for tissue in tissue_cols:
        threshold = np.percentile(merged_df[tissue], percentile)

        for substrate in substrate_types:
            # 2x2 contingency table
            is_substrate = merged_df['Predicted_Substrate'] == substrate
            is_high = merged_df[tissue] >= threshold

            # Counts
            high_substrate = (is_high & is_substrate).sum()
            high_other = (is_high & ~is_substrate).sum()
            low_substrate = (~is_high & is_substrate).sum()
            low_other = (~is_high & ~is_substrate).sum()

            # Fisher's exact test
            table = [[high_substrate, high_other],
                     [low_substrate, low_other]]

            odds_ratio, p_value = stats.fisher_exact(table)

            results.append({
                'Tissue': tissue,
                'Substrate': substrate,
                'High_Substrate': high_substrate,
                'High_Other': high_other,
                'Low_Substrate': low_substrate,
                'Low_Other': low_other,
                'Odds_Ratio': odds_ratio,
                'p_value': p_value
            })

    results_df = pd.DataFrame(results)

    # Multiple testing correction
    results_df['p_adjusted'] = results_df['p_value'] * len(results_df)
    results_df['p_adjusted'] = results_df['p_adjusted'].clip(upper=1.0)
    results_df['significant'] = results_df['p_adjusted'] < 0.05
    results_df['enriched'] = (results_df['Odds_Ratio'] > 1) & results_df['significant']
    results_df['depleted'] = (results_df['Odds_Ratio'] < 1) & results_df['significant']

    return results_df


def pairwise_substrate_comparisons(long_df: pd.DataFrame) -> pd.DataFrame:
    """Perform pairwise comparisons between substrate types per tissue."""
    logger.info("Performing pairwise Mann-Whitney U tests...")

    from itertools import combinations

    results = []
    substrates = long_df['Predicted_Substrate'].unique()

    for tissue in long_df['Tissue'].unique():
        tissue_data = long_df[long_df['Tissue'] == tissue]

        for s1, s2 in combinations(substrates, 2):
            group1 = tissue_data[tissue_data['Predicted_Substrate'] == s1]['nTPM'].values
            group2 = tissue_data[tissue_data['Predicted_Substrate'] == s2]['nTPM'].values

            if len(group1) > 0 and len(group2) > 0:
                stat, p_value = stats.mannwhitneyu(group1, group2, alternative='two-sided')

                results.append({
                    'Tissue': tissue,
                    'Comparison': f"{s1} vs {s2}",
                    'Substrate_1': s1,
                    'Substrate_2': s2,
                    'Mean_1': np.mean(group1),
                    'Mean_2': np.mean(group2),
                    'U_statistic': stat,
                    'p_value': p_value
                })

    results_df = pd.DataFrame(results)

    # Multiple testing correction
    results_df['p_adjusted'] = results_df['p_value'] * len(results_df)
    results_df['p_adjusted'] = results_df['p_adjusted'].clip(upper=1.0)
    results_df['significant'] = results_df['p_adjusted'] < 0.05

    return results_df


# =============================================================================
# 4. CO-EXPRESSION ANALYSIS
# =============================================================================

def calculate_expression_correlations(merged_df: pd.DataFrame,
                                       tissue_cols: List[str]) -> pd.DataFrame:
    """Calculate pairwise correlations of expression profiles."""
    logger.info("Calculating pairwise expression correlations...")

    # Extract expression matrix
    expr_matrix = merged_df[tissue_cols].values
    entries = merged_df['Entry'].values
    substrates = merged_df['Predicted_Substrate'].values

    n = len(entries)

    # Calculate correlation matrix (this can be memory-intensive for large datasets)
    if n > 500:
        logger.warning(f"Large dataset ({n} proteins). Correlation matrix will be {n}x{n}.")

    corr_matrix = np.corrcoef(expr_matrix)

    return corr_matrix, entries, substrates


def analyze_within_vs_between_correlation(corr_matrix: np.ndarray,
                                           substrates: np.ndarray) -> Dict:
    """Compare within-group vs between-group correlations."""
    logger.info("Analyzing within-group vs between-group expression correlations...")

    unique_substrates = np.unique(substrates)
    n = len(substrates)

    within_corrs = []
    between_corrs = []

    for i in range(n):
        for j in range(i + 1, n):
            corr = corr_matrix[i, j]
            if np.isnan(corr):
                continue

            if substrates[i] == substrates[j]:
                within_corrs.append((substrates[i], corr))
            else:
                between_corrs.append((f"{substrates[i]}-{substrates[j]}", corr))

    # Summary statistics
    within_by_type = {}
    for substrate in unique_substrates:
        type_corrs = [c for s, c in within_corrs if s == substrate]
        if type_corrs:
            within_by_type[substrate] = {
                'mean': np.mean(type_corrs),
                'median': np.median(type_corrs),
                'std': np.std(type_corrs),
                'n': len(type_corrs)
            }

    all_within = [c for _, c in within_corrs]
    all_between = [c for _, c in between_corrs]

    # Statistical test: are within-group correlations higher than between-group?
    if len(all_within) > 0 and len(all_between) > 0:
        stat, p_value = stats.mannwhitneyu(all_within, all_between, alternative='greater')
        logger.info(f"Within vs Between correlation test: U={stat:.0f}, p={p_value:.4e}")
        logger.info(f"  Mean within-group correlation: {np.mean(all_within):.4f}")
        logger.info(f"  Mean between-group correlation: {np.mean(all_between):.4f}")

    return {
        'within_by_type': within_by_type,
        'overall_within_mean': np.mean(all_within) if all_within else np.nan,
        'overall_between_mean': np.mean(all_between) if all_between else np.nan,
        'within_n': len(all_within),
        'between_n': len(all_between)
    }


# =============================================================================
# 5. VISUALIZATION HELPERS (Generate data for plotting)
# =============================================================================

def prepare_boxplot_data(long_df: pd.DataFrame, tissues: List[str] = None) -> pd.DataFrame:
    """Prepare data for box/violin plots."""
    if tissues:
        return long_df[long_df['Tissue'].isin(tissues)]
    return long_df


def prepare_tau_distribution_data(merged_df: pd.DataFrame) -> pd.DataFrame:
    """Prepare data for Tau distribution plots."""
    return merged_df[['Entry', 'Predicted_Substrate', 'Tau', 'Prediction_Confidence']].copy()


# =============================================================================
# OUTPUT GENERATION
# =============================================================================

def save_all_results(output_dir: Path, results: Dict):
    """Save all analysis results to CSV files."""
    logger.info(f"Saving results to: {output_dir}")

    # 1. Expression summary heatmap
    if 'heatmap_mean' in results:
        results['heatmap_mean'].to_csv(output_dir / 'expression_heatmap_mean.csv')
        results['heatmap_median'].to_csv(output_dir / 'expression_heatmap_median.csv')

    # 2. Expression summary statistics
    if 'expression_summary' in results:
        results['expression_summary'].to_csv(output_dir / 'expression_summary_by_substrate_tissue.csv', index=False)

    # 3. Tau analysis
    if 'tau_summary' in results:
        results['tau_summary'].to_csv(output_dir / 'tau_summary_by_substrate.csv')

    if 'merged_with_tau' in results:
        cols_to_save = ['Entry', 'Protein names', 'Predicted_Substrate',
                        'Prediction_Confidence', 'Tau', 'Max_Expression_Tissue',
                        'Max_Expression_Value']
        results['merged_with_tau'][cols_to_save].to_csv(
            output_dir / 'protein_tau_and_max_tissue.csv', index=False
        )

    if 'max_tissue_crosstab' in results:
        results['max_tissue_crosstab'].to_csv(output_dir / 'max_expression_tissue_by_substrate.csv')

    # 4. Statistical tests
    if 'kruskal_by_tissue' in results:
        results['kruskal_by_tissue'].to_csv(output_dir / 'kruskal_wallis_by_tissue.csv', index=False)

    if 'enrichment_results' in results:
        results['enrichment_results'].to_csv(output_dir / 'enrichment_high_expressors.csv', index=False)

    if 'pairwise_comparisons' in results:
        results['pairwise_comparisons'].to_csv(output_dir / 'pairwise_mannwhitney_tests.csv', index=False)

    # 5. Co-expression analysis
    if 'coexpression_summary' in results:
        coexp = results['coexpression_summary']
        coexp_df = pd.DataFrame([
            {'Metric': 'Overall within-group mean correlation', 'Value': coexp['overall_within_mean']},
            {'Metric': 'Overall between-group mean correlation', 'Value': coexp['overall_between_mean']},
            {'Metric': 'Within-group pairs count', 'Value': coexp['within_n']},
            {'Metric': 'Between-group pairs count', 'Value': coexp['between_n']},
        ])

        # Add per-substrate within-group stats
        for substrate, stats_dict in coexp['within_by_type'].items():
            coexp_df = pd.concat([coexp_df, pd.DataFrame([
                {'Metric': f'{substrate} within-group mean correlation', 'Value': stats_dict['mean']},
                {'Metric': f'{substrate} within-group pair count', 'Value': stats_dict['n']},
            ])], ignore_index=True)

        coexp_df.to_csv(output_dir / 'coexpression_summary.csv', index=False)

    # 6. Full long-format data for custom plotting
    if 'long_df' in results:
        results['long_df'].to_csv(output_dir / 'expression_long_format.csv', index=False)

    logger.info("All results saved!")


def generate_summary_report(results: Dict) -> str:
    """Generate a text summary of key findings."""
    report = []
    report.append("=" * 70)
    report.append("SAM MTase Substrate Expression Analysis - Summary Report")
    report.append("=" * 70)
    report.append("")

    # Dataset overview
    if 'merged_with_tau' in results:
        df = results['merged_with_tau']
        report.append("DATASET OVERVIEW")
        report.append("-" * 40)
        report.append(f"Total proteins analyzed: {len(df)}")
        report.append(f"Substrate type distribution:")
        for substrate, count in df['Predicted_Substrate'].value_counts().items():
            pct = count / len(df) * 100
            report.append(f"  {substrate}: {count} ({pct:.1f}%)")
        report.append("")

    # Tau analysis
    if 'tau_summary' in results:
        report.append("TISSUE SPECIFICITY (TAU) BY SUBSTRATE TYPE")
        report.append("-" * 40)
        report.append("(Tau: 0 = ubiquitous, 1 = tissue-specific)")
        report.append("")
        tau_df = results['tau_summary']
        for substrate in tau_df.index:
            row = tau_df.loc[substrate]
            report.append(f"  {substrate}:")
            report.append(f"    Mean Tau: {row['mean']:.4f} (n={int(row['count'])})")
            report.append(f"    Median Tau: {row['median']:.4f}")
        report.append("")

    # Co-expression
    if 'coexpression_summary' in results:
        coexp = results['coexpression_summary']
        report.append("CO-EXPRESSION ANALYSIS")
        report.append("-" * 40)
        report.append(f"Mean within-group correlation:  {coexp['overall_within_mean']:.4f}")
        report.append(f"Mean between-group correlation: {coexp['overall_between_mean']:.4f}")
        diff = coexp['overall_within_mean'] - coexp['overall_between_mean']
        report.append(f"Difference: {diff:.4f}")
        if diff > 0:
            report.append("→ Proteins of the same substrate type tend to have more similar expression patterns")
        report.append("")

    # Statistical tests summary
    if 'kruskal_by_tissue' in results:
        kruskal_df = results['kruskal_by_tissue']
        n_sig = kruskal_df['significant'].sum()
        report.append("STATISTICAL TESTS")
        report.append("-" * 40)
        report.append(f"Tissues with significant expression differences by substrate type:")
        report.append(f"  {n_sig} / {len(kruskal_df)} tissues (Bonferroni-adjusted p < 0.05)")
        if n_sig > 0:
            top_tissues = kruskal_df[kruskal_df['significant']].head(5)['Tissue'].tolist()
            report.append(f"  Top significant tissues: {', '.join(top_tissues)}")
        report.append("")

    if 'enrichment_results' in results:
        enrich_df = results['enrichment_results']
        enriched = enrich_df[enrich_df['enriched']]
        depleted = enrich_df[enrich_df['depleted']]
        report.append(f"Enrichment analysis (top 25% expressors per tissue):")
        report.append(f"  Significant enrichments: {len(enriched)}")
        report.append(f"  Significant depletions: {len(depleted)}")

        if len(enriched) > 0:
            report.append("  Notable enrichments (top 5 by odds ratio):")
            top_enriched = enriched.nlargest(5, 'Odds_Ratio')
            for _, row in top_enriched.iterrows():
                report.append(f"    {row['Substrate']} in {row['Tissue']}: OR={row['Odds_Ratio']:.2f}")
        report.append("")

    report.append("=" * 70)
    report.append("Analysis complete. See individual CSV files for detailed results.")
    report.append("=" * 70)

    return "\n".join(report)


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    """Main execution function."""
    logger.info("=" * 70)
    logger.info("SAM MTase Substrate Type Expression Analysis")
    logger.info("=" * 70)

    # Define paths
    script_dir = Path(__file__).parent
    expression_file = script_dir / "sam_mtase_tissue_expression_matrix.csv"

    # Look for human predictions file
    prediction_dir = script_dir.parent / "reference_files" / "prediction_files"

    # Try different possible naming conventions
    possible_names = [
        "Homo_sapiens_(Humans)_predictions.csv",
        "Homo_sapiens_predictions.csv",
        "Human_predictions.csv",
        "Homo_sapiens_(Human)_predictions.csv",
    ]

    prediction_file = None
    for name in possible_names:
        candidate = prediction_dir / name
        if candidate.exists():
            prediction_file = candidate
            break

    if prediction_file is None:
        # List available prediction files to help user
        available = list(prediction_dir.glob("*predictions.csv"))
        logger.error(f"Could not find human predictions file in: {prediction_dir}")
        logger.error(f"Looked for: {possible_names}")
        logger.error(f"Available files ({len(available)} total):")
        for f in available[:10]:
            logger.error(f"  - {f.name}")
        if len(available) > 10:
            logger.error(f"  ... and {len(available) - 10} more")
        return

    output_dir = script_dir / "substrate_expression_analysis"
    output_dir.mkdir(exist_ok=True)

    # Load data
    expression_df = load_expression_matrix(str(expression_file))
    predictions_df = load_predictions(str(prediction_file))

    # Get tissue columns (all except 'Protein names')
    tissue_cols = [c for c in expression_df.columns if c != 'Protein names']

    # Merge data
    merged_df = merge_expression_and_predictions(expression_df, predictions_df)

    if len(merged_df) == 0:
        logger.error("No proteins matched between expression and prediction data!")
        return

    # Initialize results dictionary
    results = {}

    # ==========================================================================
    # 1. Expression Profile Comparison
    # ==========================================================================
    logger.info("\n" + "=" * 50)
    logger.info("1. EXPRESSION PROFILE COMPARISON BY SUBSTRATE TYPE")
    logger.info("=" * 50)

    expression_summary, long_df = calculate_expression_summary(merged_df, tissue_cols)
    heatmap_mean, heatmap_median = create_expression_heatmap_data(expression_summary)

    results['expression_summary'] = expression_summary
    results['heatmap_mean'] = heatmap_mean
    results['heatmap_median'] = heatmap_median
    results['long_df'] = long_df

    logger.info("Expression summary calculated.")

    # ==========================================================================
    # 2. Tissue Specificity Analysis
    # ==========================================================================
    logger.info("\n" + "=" * 50)
    logger.info("2. TISSUE SPECIFICITY ANALYSIS (TAU)")
    logger.info("=" * 50)

    merged_with_tau = calculate_tau_for_all(merged_df, tissue_cols)
    tau_summary = compare_tau_by_substrate(merged_with_tau)
    max_tissue_crosstab = analyze_max_expression_tissue(merged_with_tau)

    results['merged_with_tau'] = merged_with_tau
    results['tau_summary'] = tau_summary
    results['max_tissue_crosstab'] = max_tissue_crosstab

    logger.info("\nTau summary by substrate type:")
    logger.info(f"\n{tau_summary}")

    # ==========================================================================
    # 3. Statistical Enrichment Tests
    # ==========================================================================
    logger.info("\n" + "=" * 50)
    logger.info("3. STATISTICAL ENRICHMENT TESTS")
    logger.info("=" * 50)

    kruskal_results = test_expression_differences_per_tissue(long_df)
    enrichment_results = test_enrichment_in_high_expressors(merged_with_tau, tissue_cols)
    pairwise_results = pairwise_substrate_comparisons(long_df)

    results['kruskal_by_tissue'] = kruskal_results
    results['enrichment_results'] = enrichment_results
    results['pairwise_comparisons'] = pairwise_results

    # ==========================================================================
    # 4. Co-expression Analysis
    # ==========================================================================
    logger.info("\n" + "=" * 50)
    logger.info("4. CO-EXPRESSION ANALYSIS")
    logger.info("=" * 50)

    corr_matrix, entries, substrates = calculate_expression_correlations(merged_with_tau, tissue_cols)
    coexp_summary = analyze_within_vs_between_correlation(corr_matrix, substrates)

    results['coexpression_summary'] = coexp_summary

    # ==========================================================================
    # Save Results
    # ==========================================================================
    logger.info("\n" + "=" * 50)
    logger.info("SAVING RESULTS")
    logger.info("=" * 50)

    save_all_results(output_dir, results)

    # Generate and print summary report
    report = generate_summary_report(results)
    logger.info("\n" + report)

    # Save report
    with open(output_dir / "analysis_summary_report.txt", 'w') as f:
        f.write(report)

    logger.info(f"\nAll results saved to: {output_dir}")


if __name__ == "__main__":
    main()
