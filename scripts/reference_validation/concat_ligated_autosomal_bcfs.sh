#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/slurm-%j.out

# This script concatenates the autosomal bcfs obtained from after the GLIMPSE ligation step.

# Load the bcftools module
module load samtools

# Loop trough each directory containing the ligated bcf files
for sample_cov in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/ref_panel_validation/GLIMPSE_ligate) ; do

    # Define the input directory
    INPUT_DIR=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/ref_panel_validation/GLIMPSE_ligate/${sample_cov}
    echo "Input directory of VCFs to be ligated: ${INPUT_DIR}"

    # Concatenate all the final BCF files from the autosomes (excluding the X chromosome)
    bcftools concat $(find ${INPUT_DIR}  -name "*.bcf" | grep -v "ChrX" | grep -v "autosomes") -O b -o ${INPUT_DIR}/${sample_cov}_autosomes_merged.bcf
    bcftools index ${INPUT_DIR}/${sample_cov}_autosomes_merged.bcf

done
