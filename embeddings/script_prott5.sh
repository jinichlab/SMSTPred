#!/bin/bash
#SBATCH --job-name=prot5_embeddings               # Job name
#SBATCH --ntasks=1
#SBATCH --output=embeddings_output.log       # Output file
#SBATCH --gpus=1                     # Request 1 GPU
#SBATCH --cpus-per-task=20               # Request 20 CPUs
#SBATCH --mem=64G                       # Request 50 GB of RAM, note for Prot5 more memory might be needed

#!/bin/bash

# Run your commands sequentially

# Step 1: Run the first script to create the FASTA file
python batched_prott5_multiple_files.py