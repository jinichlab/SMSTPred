import pandas as pd

df = pd.read_csv('UP_Swiss_Verified_ChEBI_SM_smiles.csv')

# Function to convert DataFrame to FASTA
def dataframe_to_fasta(df, entry_col='Entry', sequence_col='Sequence', output_file='output.fasta'):
    with open(output_file, 'w') as f:
        for _, row in df.iterrows():
            entry = row[entry_col]
            sequence = row[sequence_col]
            # Write in FASTA format
            f.write(f">{entry}\n")
            f.write(f"{sequence}\n")

# Usage
dataframe_to_fasta(df, entry_col='Entry', sequence_col='Sequence', output_file='example.fasta')
