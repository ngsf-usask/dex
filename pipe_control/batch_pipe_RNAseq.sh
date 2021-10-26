#!/bin/bash

# Runs on plato node to process RNAseq data
# Combines 4 lanes of data from nextSeq run into a single file
# Uses fastp to trim and QC reads
# Aligns to indexed genome using STAR

#SBATCH --job-name="NGSF_RNA"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0:10:00
#SBATCH --mem=4G

set -eux

THREADS=4
library=$1
raw_data_files=$2

NGSF_tag="#NGSF"

# Will echo changes when steps are completed
echo "$NGSF_tag-START True"
echo "$NGSF_tag-LIB $library"
echo "$NGSF_tag-COMBINED False"
echo "$NGSF_tag-FASTP False"
echo "$NGSF_tag-STAR False"
echo "$NGSF_tag-SBATCH False"

# activate python and then create virtual env
module load python/3.8.10
# virtual_env_path=/cluster/pipe_RNA/pipe_control/requirements.txt
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r requirements.txt # TODO MAY NEED TO ADJUST PATH IN FUTURE
echo $(python -V)


# Combine data from all four lanes
echo "$NGSF_tag -== BEGIN COMBINING FOUR LANES OF NEXTSEQ DATA ==-"
cat ${raw_data_files}${library}_*_R1_* > ${SLURM_TMPDIR}/${library}_R1.fastq.gz
cat ${raw_data_files}${library}_*_R2_* > ${SLURM_TMPDIR}/${library}_R2.fastq.gz

echo "$NGSF_tag-R1 ${library} reads: R1: $(wc -l ${SLURM_TMPDIR}/${library}_R1.fastq.gz)"
echo "$NGSF_tag-R2 ${library} reads: R2: $(wc -l ${SLURM_TMPDIR}/${library}_R2.fastq.gz)"
echo "$NGSF_tag-COMBINED True"


# Copy fastp from datastore onto node, provide permissions, and run
cp /datastore/NGSF001/software/bin/fastp ${SLURM_TMPDIR}
chmod u+x ${SLURM_TMPDIR}/fastp
ls -l


${SLURM_TMPDIR}/fastp --in1 {SLURM_TMPDIR}/${library}_R1.fastq.gz \
    --in2 {SLURM_TMPDIR}/${library}_R2.fastq.gz \
    --out1 ${library}_R1_trimmed.fastq.gz \
    --out2 ${library}_R2_trimmed.fastq.gz \
    --unpaired1 ${library}_R1_trimmed_unpaired.fastq.gz \
    --unpaired2 ${library}_R2_trimmed_unpaired.fastq.gz \
    -V \
    -l 30 \
    -p \
    -w $THREADS \
    -j ${library}.trim-report.json \
    -h ${library}.trim-report.html \
    -e 10 \
    -q 10 \
    -M 10 \
    -r \
    -W 6 \
    -g


# Move all the report documents out
# TODO create a new folder from these
mv *.json $OUTDIR
mv *.html $OUTDIR

echo "$NGSF_tag-end -== COMPLETE ==-"