#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 16
#$ -cwd

source ~/miniconda3/etc/profile.d/conda.sh
conda activate gbs

snakemake -j 16
