#!/usr/bin/env python3
"""
Run the protein-sol solubility prediction pipeline on SAM-dependent
methyltransferases and produce a results CSV + histogram.
"""

import os
import subprocess
import sys
import pandas as pd
import matplotlib.pyplot as plt

PROTEIN_SOL_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_CSV = os.path.expanduser(
    "~/Code/crossattention/files/SAM_swiss_small_molecule_smiles.csv"
)
TEMP_FASTA = "sam_mtase_input.fasta"
OUTPUT_CSV = os.path.join(PROTEIN_SOL_DIR, "solubility_results.csv")
OUTPUT_PNG = os.path.join(PROTEIN_SOL_DIR, "solubility_histogram.png")


def main():
    # 1. Read input CSV
    df = pd.read_csv(INPUT_CSV)
    print(f"Read {len(df)} entries from input CSV")

    # 2. Write temporary FASTA
    fasta_path = os.path.join(PROTEIN_SOL_DIR, TEMP_FASTA)
    with open(fasta_path, "w") as f:
        for _, row in df.iterrows():
            f.write(f">{row['Entry']}\n{row['Sequence']}\n")
    print(f"Wrote FASTA to {fasta_path}")

    # 3. Run the pipeline
    print("Running protein-sol pipeline (this may take a few minutes)...")
    result = subprocess.run(
        ["bash", "multiple_prediction_wrapper_export.sh", TEMP_FASTA],
        cwd=PROTEIN_SOL_DIR,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"Pipeline STDERR:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)
    print("Pipeline finished")

    # 4. Parse seq_prediction.txt
    pred_file = os.path.join(PROTEIN_SOL_DIR, "seq_prediction.txt")
    records = []
    with open(pred_file) as f:
        for line in f:
            if line.startswith("SEQUENCE PREDICTIONS,"):
                parts = [p.strip() for p in line.split(",")]
                # parts: [SEQUENCE PREDICTIONS, >ID, percent-sol, scaled-sol, population-sol, pI]
                entry = parts[1].lstrip(">").strip()
                scaled_sol = float(parts[3])
                records.append({"Entry": entry, "scaled_sol": scaled_sol})
    pred_df = pd.DataFrame(records)
    print(f"Parsed {len(pred_df)} predictions")

    # 5. Join with original SMILES data
    merged = pred_df.merge(df[["Entry", "smiles"]], on="Entry", how="left")
    missing = merged["smiles"].isna().sum()
    if missing > 0:
        print(f"Warning: {missing} entries missing SMILES after merge")

    # 6. Save output CSV
    merged.to_csv(OUTPUT_CSV, index=False)
    print(f"Saved results to {OUTPUT_CSV}")

    # 7. Plot histogram
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(merged["scaled_sol"], bins=30, edgecolor="black", alpha=0.7)
    ax.set_xlabel("Scaled Solubility Score")
    ax.set_ylabel("Count")
    ax.set_title("Protein-Sol Predicted Solubility (SAM MTases)")
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=150)
    print(f"Saved histogram to {OUTPUT_PNG}")

    # 8. Clean up temporary files
    for fname in [TEMP_FASTA, f"{TEMP_FASTA}_ORIGINAL"]:
        p = os.path.join(PROTEIN_SOL_DIR, fname)
        if os.path.exists(p):
            os.remove(p)
    print("Cleaned up temporary files")


if __name__ == "__main__":
    main()
