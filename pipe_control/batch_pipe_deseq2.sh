#!/bin/bash

# Runs on deseq2 using htseq output on plato, but using the rpy2 module

#SBATCH --job-name="htseq2"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00
#SBATCH --mem=20G

# conda activate diffexpr

python /globalhome/arb594/HPC/cluster/pipe_RNA/pipe_control/deseq2_htseq.py