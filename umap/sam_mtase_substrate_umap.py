#!/usr/bin/env python3
"""
SAM MTase UMAP Visualization by Predicted Substrate Type

This script creates a UMAP projection of SAM methyltransferases from the
reference proteomes (the same SAMs used to build iTOL figures), colored
by their predicted substrate binding type (RNA, DNA, Protein, Other).

Requirements:
    pip install umap-learn datamapplot pandas numpy h5py tqdm matplotlib
"""

import os
import glob
import numpy as np
import pandas as pd
import h5py
from pathlib import Path
from tqdm import tqdm
from datetime import datetime

import umap
import datamapplot
import matplotlib.pyplot as plt

# ============================================================================
# Configuration
# ============================================================================

# Paths relative to substrate_type_predictor_testing directory
SCRIPT_DIR = Path(__file__).parent.resolve()
PREDICTOR_DIR = SCRIPT_DIR.parent / "substrate_type_predictor_testing"
PRED_DIR = PREDICTOR_DIR / "reference_files" / "prediction_files"
EMBED_DIR = PREDICTOR_DIR / "reference_files" / "embeddings"
OUTPUT_DIR = SCRIPT_DIR / "sam_umap_output"

# Substrate types and their colors
SUBSTRATE_COLORS = {
    "RNA": "#cb181d",      # Red
    "DNA": "#6a51a3",      # Purple
    "Protein": "#238b45",  # Green
    "Other": "#e6550d",    # Orange
    "Mixed": "#636363",    # Gray (for multi-label cases)
}

# Random seed for reproducibility
RANDOM_SEED = 42

# ============================================================================
# Data Loading Functions
# ============================================================================

def load_embeddings_for_entries(h5_path: Path, entries: list) -> dict:
    """Load embeddings from H5 file for specific entries."""
    embeddings = {}

    if not h5_path.exists():
        return embeddings

    with h5py.File(h5_path, "r") as f:
        for entry in entries:
            if entry in f:
                emb = np.array(f[entry])
                # Validate embedding
                if emb.ndim == 1 and emb.shape[0] in [1024, 2048]:
                    if not (np.isnan(emb).any() or np.isinf(emb).any()):
                        embeddings[entry] = emb

    return embeddings


def get_predicted_substrate(row, threshold=0.5):
    """
    Determine the predicted substrate type from probability columns.
    Returns the substrate with highest probability above threshold,
    or 'Mixed' if multiple are above threshold, or the highest overall.
    """
    substrates = ["RNA", "DNA", "Protein", "Other"]
    probs = {s: row.get(s, 0) for s in substrates}

    # Convert to float, handling any non-numeric values
    for s in substrates:
        try:
            probs[s] = float(probs[s])
        except (ValueError, TypeError):
            probs[s] = 0.0

    # Find substrates above threshold
    above_threshold = [s for s, p in probs.items() if p > threshold]

    if len(above_threshold) == 1:
        return above_threshold[0]
    elif len(above_threshold) > 1:
        # Return the one with highest probability
        return max(above_threshold, key=lambda s: probs[s])
    else:
        # Nothing above threshold - return highest probability
        return max(probs, key=probs.get)


def load_all_predictions_with_embeddings():
    """
    Load all prediction files and their corresponding embeddings.
    Returns a DataFrame with entries, predictions, embeddings, and metadata.
    """
    pred_files = list(PRED_DIR.glob("*_predictions.csv"))
    print(f"Found {len(pred_files)} prediction files in {PRED_DIR}")

    all_data = []

    for pred_file in tqdm(pred_files, desc="Loading predictions and embeddings"):
        try:
            # Get organism name from filename
            organism_name = pred_file.stem.replace("_predictions", "")

            # Load predictions
            df = pd.read_csv(pred_file, low_memory=False)

            if "Entry" not in df.columns:
                continue

            # Load corresponding embeddings
            h5_path = EMBED_DIR / f"{organism_name}_embeddings.h5"
            entries = df["Entry"].tolist()
            embeddings = load_embeddings_for_entries(h5_path, entries)

            if not embeddings:
                continue

            # Process each row
            for _, row in df.iterrows():
                entry = row.get("Entry")
                if entry not in embeddings:
                    continue

                # Get predicted substrate type
                substrate = get_predicted_substrate(row)

                # Collect metadata for hover
                protein_name = row.get("Protein names", "N/A")
                organism = row.get("Organism", organism_name)
                gene_names = row.get("Gene Names", "N/A")

                # Get probabilities for hover text
                rna_prob = row.get("RNA", 0)
                dna_prob = row.get("DNA", 0)
                protein_prob = row.get("Protein", 0)
                other_prob = row.get("Other", 0)

                all_data.append({
                    "Entry": entry,
                    "Protein names": protein_name,
                    "Organism": organism,
                    "Gene Names": gene_names,
                    "Substrate": substrate,
                    "RNA_prob": rna_prob,
                    "DNA_prob": dna_prob,
                    "Protein_prob": protein_prob,
                    "Other_prob": other_prob,
                    "Embedding": embeddings[entry],
                })

        except Exception as e:
            print(f"Error processing {pred_file.name}: {e}")
            continue

    if not all_data:
        raise ValueError("No valid data found!")

    result_df = pd.DataFrame(all_data)
    print(f"Loaded {len(result_df)} SAM MTases with valid embeddings")

    return result_df


# ============================================================================
# UMAP and Visualization Functions
# ============================================================================

def compute_umap(embeddings_array, n_neighbors=15, min_dist=0.1, metric="cosine"):
    """Compute 2D UMAP projection."""
    print(f"Computing UMAP with {len(embeddings_array)} samples...")

    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        random_state=RANDOM_SEED,
        verbose=True
    )

    embedding_2d = reducer.fit_transform(embeddings_array)
    print("UMAP complete!")

    return embedding_2d


def create_hover_text(df):
    """Create hover text for each point."""
    hover_texts = []

    for _, row in df.iterrows():
        text = (
            f"Entry: {row['Entry']}\n"
            f"Protein: {row['Protein names']}\n"
            f"Organism: {row['Organism']}\n"
            f"Gene: {row['Gene Names']}\n"
            f"Predicted: {row['Substrate']}\n"
            f"---\n"
            f"RNA: {float(row['RNA_prob']):.3f}\n"
            f"DNA: {float(row['DNA_prob']):.3f}\n"
            f"Protein: {float(row['Protein_prob']):.3f}\n"
            f"Other: {float(row['Other_prob']):.3f}"
        )
        hover_texts.append(text)

    return np.array(hover_texts)


def get_custom_legend(color_mapping):
    """Generate custom CSS and HTML for the substrate type legend."""
    custom_css = """
    .row {
        display : flex;
        align-items : center;
    }
    .box {
        height:12px;
        width:12px;
        border-radius:3px;
        margin-right:8px;
    }
    #legend {
        position: absolute;
        bottom: 0;
        right: 0;
        margin: 16px;
        padding: 14px;
        border-radius: 12px;
        z-index: 2;
        background: #ffffffdd;
        font-family: Arial, sans-serif;
        font-size: 11pt;
        box-shadow: 2px 3px 10px #aaaaaa44;
    }
    #title-container {
        max-width: 75%;
    }
    """

    custom_html = '<div id="legend"><b>Substrate Type</b><br><br>\n'
    for substrate, color in color_mapping.items():
        custom_html += f'    <div class="row"><div id="{substrate}" class="box" style="background-color:{color};"></div>{substrate}</div>\n'
    custom_html += "</div>"

    return custom_css, custom_html


def create_interactive_plot(df, output_path):
    """Create and save interactive datamapplot visualization."""
    # Prepare data
    data_map = df[["UMAP_1", "UMAP_2"]].to_numpy().astype("float32")
    hover_data = create_hover_text(df)
    labels = df["Substrate"].to_numpy().astype(str)

    # Create color array based on substrate type
    color_array = df["Substrate"].map(SUBSTRATE_COLORS).fillna("#636363").values

    # Get custom legend
    custom_css, custom_html = get_custom_legend(SUBSTRATE_COLORS)

    # Create plot
    print("Creating interactive plot...")
    plot = datamapplot.create_interactive_plot(
        data_map,
        labels,
        font_family="Arial",
        title="SAM Methyltransferases by Predicted Substrate Type",
        sub_title=f"UMAP projection of {len(df):,} SAM MTases from reference proteomes",
        hover_text=hover_data,
        enable_search=True,
        marker_color_array=color_array,
        custom_css=custom_css,
        custom_html=custom_html,
        marker_size_array=np.ones(len(df)) * 5,
    )

    plot.save(str(output_path))
    print(f"Saved interactive plot to: {output_path}")

    return plot


def create_static_plot(df, output_path):
    """Create and save a static matplotlib visualization."""
    fig, ax = plt.subplots(figsize=(14, 10))

    # Plot each substrate type
    for substrate, color in SUBSTRATE_COLORS.items():
        mask = df["Substrate"] == substrate
        if mask.sum() > 0:
            ax.scatter(
                df.loc[mask, "UMAP_1"],
                df.loc[mask, "UMAP_2"],
                c=color,
                label=f"{substrate} (n={mask.sum():,})",
                alpha=0.6,
                s=8,
                edgecolors="none"
            )

    ax.set_xlabel("UMAP 1", fontsize=12)
    ax.set_ylabel("UMAP 2", fontsize=12)
    ax.set_title("SAM Methyltransferases by Predicted Substrate Type", fontsize=14, fontweight="bold")
    ax.legend(loc="upper right", fontsize=10, framealpha=0.9)

    # Remove top and right spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Saved static plot to: {output_path}")


# ============================================================================
# Main
# ============================================================================

def main():
    print("=" * 70)
    print("SAM MTase UMAP Visualization by Predicted Substrate Type")
    print("=" * 70)

    # Create output directory
    OUTPUT_DIR.mkdir(exist_ok=True)

    # Load data
    print("\n[1/4] Loading predictions and embeddings...")
    df = load_all_predictions_with_embeddings()

    # Print substrate distribution
    print("\nSubstrate type distribution:")
    for substrate in SUBSTRATE_COLORS.keys():
        count = (df["Substrate"] == substrate).sum()
        pct = count / len(df) * 100
        print(f"  {substrate}: {count:,} ({pct:.1f}%)")

    # Stack embeddings
    print("\n[2/4] Preparing embeddings...")
    embeddings_array = np.vstack(df["Embedding"].values).astype(np.float32)
    print(f"Embeddings shape: {embeddings_array.shape}")

    # Compute UMAP
    print("\n[3/4] Computing UMAP projection...")
    umap_coords = compute_umap(embeddings_array)
    df["UMAP_1"] = umap_coords[:, 0]
    df["UMAP_2"] = umap_coords[:, 1]

    # Save DataFrame (without embeddings for size)
    date_str = datetime.now().strftime("%Y-%m-%d")
    csv_path = OUTPUT_DIR / f"sam_mtase_umap_{date_str}.csv"
    df_save = df.drop(columns=["Embedding"])
    df_save.to_csv(csv_path, index=False)
    print(f"\nSaved DataFrame to: {csv_path}")

    # Create visualizations
    print("\n[4/4] Creating visualizations...")

    # Static plot
    static_path = OUTPUT_DIR / f"sam_mtase_umap_static_{date_str}.png"
    create_static_plot(df, static_path)

    # Interactive plot
    interactive_path = OUTPUT_DIR / f"sam_mtase_umap_interactive_{date_str}.html"
    create_interactive_plot(df, interactive_path)

    # Summary
    print("\n" + "=" * 70)
    print("COMPLETE!")
    print("=" * 70)
    print(f"\nOutput files in: {OUTPUT_DIR}")
    print(f"  - {csv_path.name} (data)")
    print(f"  - {static_path.name} (static image)")
    print(f"  - {interactive_path.name} (interactive HTML)")
    print("\nOpen the HTML file in a browser for interactive exploration!")


if __name__ == "__main__":
    main()
