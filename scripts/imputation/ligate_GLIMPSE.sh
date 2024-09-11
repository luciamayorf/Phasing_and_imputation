#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/slurm-%j.out

# This script performs the ligation of the imputed bcfs (separated by chunks) with GLIMPSE (Third step). Usage of the script: 
# for chr_dir in $(ls -d /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_all/GLIMPSE_phase/*); do
#   sbatch -t 01:00:00 --mem 3GB ligate_GLIMPSE.sh ${chr_dir} <imputation_batchID>
# done

# Load the GLIMPSE1.1 module
module load cesga/2020 gcc/system glimpse/1.1.1

# Define the input directory
INPUT_DIR=${1}
echo "Input directory of ligating step: ${INPUT_DIR}"

# Define imputation batch of samples ID
ID_IMP=${2}
echo "Batch imputation ID: ${ID_IMP}"

# Define the chromosome
CHR=$(echo ${INPUT_DIR} | awk -F'/' '{print $NF}')
echo "Chromosome to be phased: ${CHR}"

# Define the output directory
OUT_DIR=$(echo "${INPUT_DIR}" | sed "s/GLIMPSE_phase\/${CHR}/GLIMPSE_ligate/g")
mkdir -p ${OUT_DIR}
echo "Output directory of ligating step: ${OUT_DIR}"

# Define list of BCF files to be ligated
find ${INPUT_DIR} -name "*.bcf" > ${INPUT_DIR}/${CHR}_bcf_list.txt

# Ligate the BCF files and index the output
GLIMPSE_ligate --input ${INPUT_DIR}/${CHR}_bcf_list.txt --output ${OUT_DIR}/${ID_IMP}_${CHR}_merged.bcf 
bcftools index -f ${OUT_DIR}/${ID_IMP}_${CHR}_merged.bcf 

# Delete the list of BCF files
rm ${INPUT_DIR}/${CHR}_bcf_list.txt

