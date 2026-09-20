#!/bin/bash

#PBS -l select=1:ncpus=8:ngpus=1:mem=32gb
#PBS -l walltime=3:00:00
#PBS -m ae
#PBS -j oe
#PBS -N colabfold

# Usage:
# qsub -v INPUT_FASTA=motA_top10_distant.fa,OUTPUT_DIR=output_motA ColabFold_batch_PBS.sh

# Move to the directory where the job was submitted
cd "$PBS_O_WORKDIR" || exit

# Activate Conda environment
conda init
conda activate myenv

# Run ColabFold
colabfold_batch --num-recycle 5 "$INPUT_FASTA" "$OUTPUT_DIR"
