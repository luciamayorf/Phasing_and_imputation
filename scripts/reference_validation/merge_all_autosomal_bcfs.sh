#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/slurm-%j.out

# This script merges the autosomal bcfs obtained from the whole imputation process.

# Load the bcftools module
module load samtools

# Define the input directory
INPUT_DIR=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/ref_panel_validation/GLIMPSE_ligate

# Loop trough each coverage level to merge all the obtained BCF files
for cov in 0_5 1 2; do

    # Merge all the final BCF files from the autosomes (excluding the X chromosome)
    bcftools merge $(find ${INPUT_DIR}  -name "*_${cov}x_autosomes_merged.bcf") -O b -o ${INPUT_DIR}/c_lp_all_samples_${cov}x_autosomes_imputed.bcf

done

