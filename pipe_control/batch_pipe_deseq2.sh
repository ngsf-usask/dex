#!/bin/bash

# Installs deseq2 on the node, and runs using python script
# https://github.com/wckdouglas/diffexpr

#SBATCH --job-name="NGSF_RNA"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0:20:00
#SBATCH --mem=4G

set -eux

# conda activate diffexpr


# Transfer files, and verify installation
mkdir ${SLURM_TMPDIR}/diffexpr
rsync -rvzP /datastore/NGSF001/software/src/diffexpr/diffexpr-master/* ${SLURM_TMPDIR}/diffexpr/


# Use direct paths to env created in conda 
# /globalhome/arb594/HPC/anaconda3/envs/diffexpr/bin/Rscript ${SLURM_TMPDIR}/diffexpr/setup.R
# /globalhome/arb594/HPC/anaconda3/envs/diffexpr/bin/python ${SLURM_TMPDIR}/diffexpr/setup.py install


/globalhome/arb594/HPC/anaconda3/envs/diffexpr/bin/python /globalhome/arb594/HPC/cluster/pipe_RNA/pipe_control/deseq2_htseq.py
# TODO adjust path in future