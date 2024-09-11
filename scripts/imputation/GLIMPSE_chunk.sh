#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/slurm-%j.out

# This script generates the chunks for the GLIMPSE imputation (First step). It only needs to be run once (if the reference panel doesn't change). 
# The idea is to launch one job per chromosome. Usage of the script: 
# for chr in $(cut -f1 /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.big_chromosomes.bed); do
#   sbatch  -t 00:30:00 --mem 3GB chunks_GLIMPSE.sh ${chr} <input_VCFref>
# done

# Load the GLIMPSE1.1 module
module load cesga/2020 gcc/system glimpse/1.1.1

# Define the chromosome
CHR=${1}

# Define the input VCF reference file
VCF_REF=${2}

# Define the output directory
OUT_DIR=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/GLIMPSE_chunk
mkdir -p ${OUT_DIR}

# Define output name
OUT_NAME=$(basename ${VCF_REF} .vcf.gz)

# Run GLIMPSE
GLIMPSE_chunk --input ${VCF_REF} --region ${CHR} --window-size 2000000 --buffer-size 200000 --output ${OUT_DIR}/${OUT_NAME}_${CHR}_chunks.txt
