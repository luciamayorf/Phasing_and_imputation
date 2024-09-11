#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/slurm-%j.out

# This script performs the sampling (defining haplotypes) of the imputed and ligated bcfs with GLIMPSE (Fourth step). Usage of the script: 
# for bcf in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_all/GLIMPSE_ligate/*.bcf); do
#   sbatch -t 00:30:00 --mem 2GB sample_GLIMPSE.sh ${bcf}
# done

# Load the GLIMPSE1.1 module
module load cesga/2020 gcc/system glimpse/1.1.1

# Define the input BCF file
INPUT_BCF=${1}
echo "Input BCF file of sampling step: ${INPUT_BCF}"

# Define the output directory
OUT_DIR=$(dirname ${INPUT_BCF} | sed "s/GLIMPSE_ligate/GLIMPSE_sample/g")
mkdir -p ${OUT_DIR}
echo "Output directory of sampling step: ${OUT_DIR}"

# Define the output VCF name
OUT_BCF=$(basename ${INPUT_BCF} | sed 's/_merged.bcf/_phased.bcf/g')
echo "Output BCF file of sampling step: ${OUT_BCF}"

# Sample haplotypes from the chr BCF files and index the output
GLIMPSE_sample --input ${INPUT_BCF} --solve --output ${OUT_DIR}/${OUT_BCF} 
bcftools index -f ${OUT_DIR}/${OUT_BCF} 
