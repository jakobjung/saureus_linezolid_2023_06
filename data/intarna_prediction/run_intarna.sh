#!/bin/bash
#SBATCH --job-name=intarna_job
#SBATCH --output=intarna_output.txt
#SBATCH --error=intarna_error.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G

# Load any necessary modules or environments
# module load ...

# Run the IntaRNA command
IntaRNA -q sig_srnas.fasta -t all_transcripts.fasta --outMode=C --out intarna_result.csv -v
