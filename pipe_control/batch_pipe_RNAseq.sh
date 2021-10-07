#!/bin/bash

# Runs on plato node to process RNAseq data
# Combines 4 lanes of data from nextSeq run into a single file
# Uses fastp to trim and QC reads
# Aligns to indexed genome using STAR

#SBATCH --job-name=pipe_RNA
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0:5:00
#SBATCH --mem=4G

set -eux

THREADS=4
library=$1
raw_data_files=$2

NGSF_tag="#NGSF"

echo "$NGSF_tag JOBID $SLURM_JOBID"

# activate python and then create virtual env
module load python/3.8.10
# virtual_env_path=/cluster/pipe_RNA/pipe_control/requirements.txt
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r requirements.txt # TODO MAY NEED TO ADJUST PATH IN FUTURE
echo $(python -V)


# Combine data from all four lanes
echo "-== BEGIN COMBINING FOUR LANES OF NEXTSEQ DATA ==-"
cat ${raw_data_files}${library}_*_R1_* > ${SLURM_TMPDIR}/${library}_R1.fastq.gz
cat ${raw_data_files}${library}_*_R2_* > ${SLURM_TMPDIR}/${library}_R2.fastq.gz
echo "$NGSF_tag ${library} reads: R1: $(wc -l ${SLURM_TMPDIR}/${library}_R1.fastq.gz)"
echo "$NGSF_tag ${library} reads: R2: $(wc -l ${SLURM_TMPDIR}/${library}_R2.fastq.gz)"




echo "$NGSF_tag-end -== COMPLETE ==-"