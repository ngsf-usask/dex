#!/bin/bash

# Runs on plato node to process RNAseq data
# Combines 4 lanes of data from nextSeq run into a single file
# Uses fastp to trim and QC reads
# Aligns to indexed genome using STAR

#SBATCH --job-name="NGSF_RNA"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=40G

set -eux

THREADS=4
library=$1
raw_data_files=$2
index=$3
outdir=$4
gtf=$5

NGSF_tag="#NGSF"

# Will echo changes when steps are completed
echo "$NGSF_tag-RUN True"
echo "$NGSF_tag-LIB $library"
echo "$NGSF_tag-COMBINED False"
echo "$NGSF_tag-FASTP False"
echo "$NGSF_tag-STAR False"
echo "$NGSF_tag-HTSEQ False"
echo "$NGSF_tag-SBATCH False"

# activate python and then create virtual env
module load python/3.8.10
# virtual_env_path=/cluster/pipe_RNA/pipe_control/requirements.txt
virtualenv --no-download ${SLURM_TMPDIR}/env
source ${SLURM_TMPDIR}/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r ~/cluster/pipe_RNA/pipe_control/requirements_htseq.txt # TODO MAY NEED TO ADJUST PATH IN FUTURE
echo $(python -V)

# Combine data from all four lanes

echo "$NGSF_tag -== BEGIN COMBINING FOUR LANES OF NEXTSEQ DATA ==-"
cat ${raw_data_files}${library}_*_R1_* > ${SLURM_TMPDIR}/${library}_R1.fastq.gz
cat ${raw_data_files}${library}_*_R2_* > ${SLURM_TMPDIR}/${library}_R2.fastq.gz

echo "$NGSF_tag-R1 ${library} reads: R1: $(wc -l ${SLURM_TMPDIR}/${library}_R1.fastq.gz)"
echo "$NGSF_tag-R2 ${library} reads: R2: $(wc -l ${SLURM_TMPDIR}/${library}_R2.fastq.gz)"

echo "$NGSF_tag-COMBINED True"


# Activate fastp and run 
module load fastp/0.20.1

fastp --in1 ${SLURM_TMPDIR}/${library}_R1.fastq.gz \
    --in2 ${SLURM_TMPDIR}/${library}_R2.fastq.gz \
    --out1 ${SLURM_TMPDIR}/${library}_R1_trimmed.fastq.gz \
    --out2 ${SLURM_TMPDIR}/${library}_R2_trimmed.fastq.gz \
    --unpaired1 ${SLURM_TMPDIR}/${library}_R1_trimmed_unpaired.fastq.gz \
    --unpaired2 ${SLURM_TMPDIR}/${library}_R2_trimmed_unpaired.fastq.gz \
    -V \
    -l 30 \
    -p \
    -w $THREADS \
    -j ${library}.trim-report.fastp.json \
    -h ${library}.trim-report.fastp.html \
    -e 10 \
    -q 10 \
    -M 10 \
    -r \
    -W 6 \
    -g

echo "$NGSF_tag-FASTP True"

# Create directory to store library output
mkdir ${SLURM_TMPDIR}/${library}

mv *.json ${SLURM_TMPDIR}/${library}
mv *.html ${SLURM_TMPDIR}/${library}

# Begin STAR analysis
module load cellranger/2.1.0
module load samtools/1.10

PATH=/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/cellranger/2.1.0/STAR/5dda596/:$PATH

# Move into library output so that output is automatically generated here
cd ${SLURM_TMPDIR}/${library}

STAR --runMode alignReads \
    --runThreadN $THREADS \
    --genomeDir $index \
    --readFilesIn <(zcat ../${library}_R1_trimmed.fastq.gz) <(zcat ../${library}_R2_trimmed.fastq.gz) \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outSAMtype BAM SortedByCoordinate

samtools index Aligned.sortedByCoord.out.bam

# rsync -rvz ${SLURM_TMPDIR}/${library} $outdir

echo "$NGSF_tag-STAR True"

htseq-count -f bam \
    -r pos \
    -s reverse \
    -m union \
    --nonunique all \
    Aligned.sortedByCoord.out.bam \
    $gtf > ${library}.counts.htseq.txt


echo "$NGSF_tag-HTSEQ True"

rsync -rvz ${SLURM_TMPDIR}/${library} $outdir

echo "$NGSF_tag-end -== COMPLETE ==-"