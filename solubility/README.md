# Solubility

Python wrapper for running [protein-sol](https://protein-sol.manchester.ac.uk/) solubility predictions on SAM-dependent methyltransferases.

## Scripts

| Script | Description |
|--------|-------------|
| `run_solubility.py` | Wrapper that converts CSV input to FASTA, runs the protein-sol prediction pipeline, parses results, and merges with SMILES data. Outputs solubility scores and a histogram. |

## Dependencies

This module requires the **protein-sol** prediction pipeline to be installed separately. The pipeline includes Perl scripts and reference data not included in this repository.

- **protein-sol**: [https://protein-sol.manchester.ac.uk/](https://protein-sol.manchester.ac.uk/)
- Hebditch et al. (2017) "protein-sol: a web tool for predicting protein solubility from sequence." *Bioinformatics*, 33(19), 3098–3100.

## Usage

```bash
# Run solubility predictions
python run_solubility.py
```

## Input/Output

- **Input**: CSV with Entry and Sequence columns
- **Output**: CSV with scaled solubility scores + SMILES, histogram PNG
