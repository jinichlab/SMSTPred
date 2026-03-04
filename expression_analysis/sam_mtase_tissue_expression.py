#!/usr/bin/env python3
"""
SAM MTase Tissue Expression Analysis

This script:
1. Filters human proteome for SAM MTases (IPR029063)
2. Maps UniProt IDs to Ensembl Gene IDs
3. Retrieves tissue expression data from Human Protein Atlas
4. Outputs a merged dataset with expression information
"""

import pandas as pd
import logging
import sys
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Constants
INTERPRO_SAM_MTASE = "IPR029063"
HPA_TISSUE_URL = "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip"


def load_proteome(filepath: str) -> pd.DataFrame:
    """Load the human proteome CSV file."""
    logger.info(f"Loading proteome from: {filepath}")
    df = pd.read_csv(filepath)
    logger.info(f"Total proteins in proteome: {len(df):,}")
    return df


def filter_sam_mtases(df: pd.DataFrame) -> pd.DataFrame:
    """Filter for SAM MTases by searching for IPR029063 in InterPro column."""
    logger.info(f"Filtering for SAM MTases (InterPro: {INTERPRO_SAM_MTASE})")

    # Handle NaN values and filter
    mask = df['InterPro'].fillna('').str.contains(INTERPRO_SAM_MTASE, regex=False)
    sam_mtases = df[mask].copy()

    logger.info(f"SAM MTases found: {len(sam_mtases)}")

    # Log some details about the proteins found
    if len(sam_mtases) > 0:
        logger.info("Sample SAM MTase proteins:")
        for _, row in sam_mtases.head(5).iterrows():
            name = row['Protein names'][:60] + '...' if len(str(row['Protein names'])) > 60 else row['Protein names']
            logger.info(f"  - {row['Entry']}: {name}")

    return sam_mtases


def map_uniprot_to_ensembl(uniprot_ids: list) -> pd.DataFrame:
    """Map UniProt IDs to Ensembl Gene IDs using UniProt ID Mapping API."""
    logger.info(f"Mapping {len(uniprot_ids)} UniProt IDs to Ensembl Gene IDs...")

    try:
        from UniProtMapper import ProtMapper
    except ImportError:
        logger.error("UniProtMapper not installed. Install with: pip install uniprot-id-mapper")
        logger.info("Attempting to install uniprot-id-mapper...")
        import subprocess
        subprocess.check_call([sys.executable, "-m", "pip", "install", "uniprot-id-mapper"])
        from UniProtMapper import ProtMapper

    mapper = ProtMapper()

    # Perform the mapping
    result, failed = mapper.get(
        ids=uniprot_ids,
        from_db="UniProtKB_AC-ID",
        to_db="Ensembl"
    )

    # Log results
    successful_mappings = len(result) if result is not None and not result.empty else 0
    failed_mappings = len(failed) if failed else 0

    logger.info(f"Successfully mapped: {successful_mappings}")
    logger.info(f"Failed to map: {failed_mappings}")

    if failed and len(failed) > 0:
        logger.warning(f"Failed IDs (first 10): {failed[:10]}")

    return result, failed


def fetch_gene_names_for_failed_ids(failed_ids: list) -> dict:
    """Fetch gene names from UniProt for IDs that failed Ensembl mapping.

    Returns:
        dict: Mapping of UniProt ID -> primary gene name
    """
    if not failed_ids:
        return {}

    logger.info(f"Fetching gene names for {len(failed_ids)} failed IDs via UniProt API...")

    import urllib.request
    import json

    gene_name_map = {}

    for uniprot_id in failed_ids:
        try:
            url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
            with urllib.request.urlopen(url, timeout=10) as response:
                data = json.loads(response.read().decode())

            # Extract primary gene name
            genes = data.get('genes', [])
            if genes:
                primary_gene = genes[0].get('geneName', {}).get('value')
                if primary_gene:
                    gene_name_map[uniprot_id] = primary_gene
                    logger.info(f"  {uniprot_id} -> gene name: {primary_gene}")
                else:
                    logger.info(f"  {uniprot_id} -> no gene name found")
            else:
                logger.info(f"  {uniprot_id} -> no gene info")

        except Exception as e:
            logger.warning(f"  {uniprot_id} -> failed to fetch: {e}")

    return gene_name_map


def download_hpa_tissue_data(output_dir: str) -> pd.DataFrame:
    """Download and load HPA tissue expression data."""
    import urllib.request
    import zipfile
    import os

    output_path = Path(output_dir)
    tsv_file = output_path / "rna_tissue_consensus.tsv"
    zip_file = output_path / "rna_tissue_consensus.tsv.zip"

    # Check if already downloaded
    if tsv_file.exists():
        logger.info(f"HPA tissue data already exists at: {tsv_file}")
    else:
        logger.info(f"Downloading HPA tissue expression data from: {HPA_TISSUE_URL}")

        # Download
        urllib.request.urlretrieve(HPA_TISSUE_URL, zip_file)
        logger.info("Download complete. Extracting...")

        # Extract
        with zipfile.ZipFile(zip_file, 'r') as z:
            z.extractall(output_path)

        # Clean up zip file
        os.remove(zip_file)
        logger.info("Extraction complete.")

    # Load the data
    logger.info("Loading HPA tissue expression data...")
    hpa_df = pd.read_csv(tsv_file, sep='\t')

    unique_genes = hpa_df['Gene'].nunique() if 'Gene' in hpa_df.columns else 0
    unique_tissues = hpa_df['Tissue'].nunique() if 'Tissue' in hpa_df.columns else 0

    logger.info(f"HPA data loaded: {len(hpa_df):,} rows")
    logger.info(f"Unique genes in HPA: {unique_genes:,}")
    logger.info(f"Unique tissues in HPA: {unique_tissues}")

    return hpa_df


def merge_expression_data(sam_mtases: pd.DataFrame,
                          mapping_result: pd.DataFrame,
                          hpa_data: pd.DataFrame,
                          failed_ids: list = None,
                          failed_gene_names: dict = None) -> tuple:
    """Merge SAM MTase data with expression data.

    Args:
        sam_mtases: DataFrame of SAM MTase proteins
        mapping_result: DataFrame with UniProt to Ensembl mappings
        hpa_data: DataFrame with HPA tissue expression data
        failed_ids: List of UniProt IDs that failed Ensembl mapping
        failed_gene_names: Dict mapping failed UniProt IDs to gene names (for fallback)

    Returns:
        tuple: (merged_data, unmatched_df) where unmatched_df contains SAM MTases
               that could not be matched to tissue expression data.
    """
    logger.info("Merging expression data...")

    # First, merge SAM MTases with their Ensembl IDs
    # The mapping result has 'From' (UniProt) and 'To' (Ensembl) columns
    if mapping_result is None or mapping_result.empty:
        logger.warning("No mapping results to merge")
        return pd.DataFrame(), sam_mtases.copy()

    # Rename columns for clarity
    mapping_df = mapping_result.rename(columns={'From': 'Entry', 'To': 'Ensembl_ID'})

    # Merge SAM MTases with Ensembl IDs
    sam_with_ensembl = sam_mtases.merge(mapping_df[['Entry', 'Ensembl_ID']], on='Entry', how='left')

    mapped_count = sam_with_ensembl['Ensembl_ID'].notna().sum()
    logger.info(f"SAM MTases with Ensembl IDs: {mapped_count}")

    # HPA uses 'Gene' column for Ensembl IDs (format: ENSG...)
    # We need to extract the Ensembl gene ID from the full ID if needed
    sam_with_ensembl['Ensembl_Gene'] = sam_with_ensembl['Ensembl_ID'].str.extract(r'(ENSG\d+)', expand=False)

    # Merge with HPA data by Ensembl Gene ID
    # HPA 'Gene' column contains Ensembl IDs
    merged = sam_with_ensembl.merge(
        hpa_data,
        left_on='Ensembl_Gene',
        right_on='Gene',
        how='left'
    )

    # Count how many SAM MTases have expression data
    with_expression = merged.dropna(subset=['nTPM'])
    unique_with_expression = with_expression['Entry'].nunique()

    logger.info(f"SAM MTases matched to tissue data (via Ensembl): {unique_with_expression}")

    # Fallback: Try to match failed IDs by gene name
    if failed_gene_names:
        logger.info("Attempting gene name fallback for failed Ensembl mappings...")

        # Get unique gene names in HPA
        hpa_gene_names = set(hpa_data['Gene name'].unique())

        recovered_count = 0
        for uniprot_id, gene_name in failed_gene_names.items():
            if gene_name in hpa_gene_names:
                # This gene name exists in HPA - we can recover it!
                logger.info(f"  RECOVERED: {uniprot_id} via gene name '{gene_name}'")

                # Get HPA data for this gene name
                hpa_subset = hpa_data[hpa_data['Gene name'] == gene_name].copy()

                # Get the SAM MTase row for this entry
                sam_row = sam_mtases[sam_mtases['Entry'] == uniprot_id].copy()

                if not sam_row.empty and not hpa_subset.empty:
                    # Create merged rows for this protein
                    sam_row['Ensembl_ID'] = None
                    sam_row['Ensembl_Gene'] = hpa_subset['Gene'].iloc[0]  # Get Ensembl ID from HPA

                    for _, hpa_row in hpa_subset.iterrows():
                        new_row = sam_row.iloc[0].to_dict()
                        new_row.update(hpa_row.to_dict())
                        merged = pd.concat([merged, pd.DataFrame([new_row])], ignore_index=True)

                    recovered_count += 1
            else:
                logger.info(f"  NOT FOUND: {uniprot_id} gene name '{gene_name}' not in HPA")

        if recovered_count > 0:
            logger.info(f"Recovered {recovered_count} proteins via gene name fallback")

            # Recount
            with_expression = merged.dropna(subset=['nTPM'])
            unique_with_expression = with_expression['Entry'].nunique()

    logger.info(f"SAM MTases matched to tissue data (total): {unique_with_expression}")
    logger.info(f"Total expression records: {len(with_expression):,}")

    # Identify unmatched SAM MTases
    matched_entries = set(with_expression['Entry'].unique())
    all_entries = set(sam_mtases['Entry'].unique())
    unmatched_entries = all_entries - matched_entries

    # Create dataframe of unmatched SAM MTases with reason
    unmatched_df = sam_mtases[sam_mtases['Entry'].isin(unmatched_entries)].copy()

    # Add reason for non-match
    def get_unmatched_reason(entry):
        if failed_ids and entry in failed_ids:
            return "No Ensembl ID mapping"
        # Check if it had an Ensembl ID but no HPA match
        entry_ensembl = sam_with_ensembl[sam_with_ensembl['Entry'] == entry]['Ensembl_Gene'].dropna()
        if len(entry_ensembl) == 0:
            return "No Ensembl ID mapping"
        else:
            return f"Ensembl ID not in HPA database ({entry_ensembl.iloc[0]})"

    unmatched_df['Unmatched_Reason'] = unmatched_df['Entry'].apply(get_unmatched_reason)

    logger.info(f"SAM MTases NOT matched to tissue data: {len(unmatched_df)}")

    return merged, unmatched_df


def create_summary_table(merged_data: pd.DataFrame, unmatched_df: pd.DataFrame, output_path: str):
    """Create a summary table of SAM MTase expression across tissues."""
    if merged_data.empty:
        logger.warning("No data to summarize")
        return

    # Create a pivot table: proteins vs tissues
    expression_data = merged_data.dropna(subset=['nTPM'])

    if expression_data.empty:
        logger.warning("No expression data found")
        return

    # Pivot to get expression matrix
    pivot_table = expression_data.pivot_table(
        index='Entry',
        columns='Tissue',
        values='nTPM',
        aggfunc='first'
    )

    # Add protein names
    protein_names = merged_data[['Entry', 'Protein names']].drop_duplicates()
    pivot_table = pivot_table.merge(protein_names, left_index=True, right_on='Entry')
    pivot_table = pivot_table.set_index('Entry')

    # Reorder columns to put protein names first
    cols = ['Protein names'] + [c for c in pivot_table.columns if c != 'Protein names']
    pivot_table = pivot_table[cols]

    # Save
    output_file = Path(output_path) / "sam_mtase_tissue_expression_matrix.csv"
    pivot_table.to_csv(output_file)
    logger.info(f"Expression matrix saved to: {output_file}")

    # Also save the full merged data
    full_output = Path(output_path) / "sam_mtase_tissue_expression_full.csv"
    merged_data.to_csv(full_output, index=False)
    logger.info(f"Full data saved to: {full_output}")

    # Save unmatched SAM MTases
    if not unmatched_df.empty:
        unmatched_output = Path(output_path) / "sam_mtase_unmatched_to_tissue_data.csv"
        # Select relevant columns for the output
        output_cols = ['Entry', 'Protein names', 'Gene Names', 'Unmatched_Reason']
        available_cols = [c for c in output_cols if c in unmatched_df.columns]
        unmatched_df[available_cols].to_csv(unmatched_output, index=False)
        logger.info(f"Unmatched SAM MTases saved to: {unmatched_output}")

        # Log the unmatched proteins
        logger.info("=" * 60)
        logger.info("UNMATCHED SAM MTases (no tissue expression data):")
        logger.info("=" * 60)
        for _, row in unmatched_df.iterrows():
            protein_name = str(row['Protein names'])[:50] + '...' if len(str(row['Protein names'])) > 50 else row['Protein names']
            logger.info(f"  {row['Entry']}: {protein_name}")
            logger.info(f"    Reason: {row['Unmatched_Reason']}")

    return pivot_table


def main():
    """Main execution function."""
    logger.info("=" * 60)
    logger.info("SAM MTase Tissue Expression Analysis")
    logger.info("=" * 60)

    # Define paths
    script_dir = Path(__file__).parent
    # Proteome file is in the reference_proteomes directory
    proteome_file = script_dir.parent / "reference_files" / "reference_proteomes" / "Homo_sapiens_(Humans).csv"
    output_dir = script_dir

    # Step 1: Load proteome
    proteome = load_proteome(str(proteome_file))

    # Step 2: Filter for SAM MTases
    sam_mtases = filter_sam_mtases(proteome)

    if len(sam_mtases) == 0:
        logger.error("No SAM MTases found. Exiting.")
        return

    # Step 3: Map UniProt to Ensembl
    uniprot_ids = sam_mtases['Entry'].tolist()
    mapping_result, failed_ids = map_uniprot_to_ensembl(uniprot_ids)

    # Step 3b: Fetch gene names for failed IDs (for fallback matching)
    failed_gene_names = {}
    if failed_ids:
        failed_gene_names = fetch_gene_names_for_failed_ids(failed_ids)

    # Step 4: Download HPA tissue data
    hpa_data = download_hpa_tissue_data(str(output_dir))

    # Step 5: Merge data (with gene name fallback for failed Ensembl mappings)
    merged_data, unmatched_df = merge_expression_data(
        sam_mtases, mapping_result, hpa_data, failed_ids, failed_gene_names
    )

    # Step 6: Create summary
    summary = create_summary_table(merged_data, unmatched_df, str(output_dir))

    # Final summary
    logger.info("=" * 60)
    logger.info("SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total proteins in proteome:     {len(proteome):,}")
    logger.info(f"SAM MTases identified:          {len(sam_mtases)}")
    logger.info(f"Successfully mapped to Ensembl: {len(mapping_result) if mapping_result is not None else 0}")
    logger.info(f"Failed to map:                  {len(failed_ids) if failed_ids else 0}")

    if not merged_data.empty:
        unique_with_expr = merged_data.dropna(subset=['nTPM'])['Entry'].nunique()
        logger.info(f"Matched to tissue data:         {unique_with_expr}")
        logger.info(f"NOT matched to tissue data:     {len(unmatched_df)}")

    logger.info("=" * 60)
    logger.info("Analysis complete!")


if __name__ == "__main__":
    main()
