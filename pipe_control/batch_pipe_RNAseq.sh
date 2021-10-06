#!/bin/bash

# Runs on plato node to process RNAseq data
# Combines 4 lanes of data from nextSeq run into a single file
# Uses fastp to trim and QC reads
# Aligns to indexed genome using STAR

#SBATCH --job-name=RNA_trim
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0:5:00
#SBATCH --mem=4G

set -eux

THREADS=4
infile=$1
raw_data_files=$2

# activate virtualenv and python

# add for loop

echo $infile
echo $raw_data_files