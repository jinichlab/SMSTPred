#!/bin/bash
#SBATCH --job-name=cosine_classifier_testing  # Job name
#SBATCH --ntasks=1
#SBATCH --output=testing_output.log           # Output file
#SBATCH --cpus-per-task=20                     # Request 20 CPUs
#SBATCH --mem=20G                              # Request 20 GB of RAM

# Load any required modules (if needed)
# module load python  # Uncomment and modify if your environment requires loading specific modules

# Run the desired command
python SAM_mlp_xgboost_comparison_mtb_predictions.py