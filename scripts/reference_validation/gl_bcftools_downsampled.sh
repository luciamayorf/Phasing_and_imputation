#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/downsampling/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/downsampling/slurm-%j.out

# In this script, we compute the Genotypes Likelihood of a BAM using bcftools mpileup and call. 
# It requires a VCF and TSV files of the reference panel to specify the target variants, as well as a reference genome to align the BAM.

# This script was specifically designed to compute the Genotypes Likelihood of the leave-one-out approach to estimate the imputation accuracy.
# The SAMPLE, VCF, TSV variable definition should be simpler for straight imputation.

# Usage of the script: sbatch gl_bcftools_downsampled.sh <input_bam> <vcfs_directory>

# Load the samtools module:
module load cesga/2020 samtools/1.19

# Define input vcf
BAM=${1}

# Define the reference genome
REFGEN=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.scaffolds.revcomp.scaffolds.fa

# Define the sample name
SAMPLE=$(basename ${BAM} | sed 's/mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.*x.bam//')
echo "Sample BAM ID: ${SAMPLE}"

# Define coverage of the downsampled BAM:
COV=$(basename ${BAM} .bam | sed "s/${SAMPLE}mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.//")
echo "Coverage: ${COV}"

# Define the new sample name
SAMPLE_ID=$(grep "${SAMPLE}" /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/FASTQ_files/novogene_lp_sept2023/fastq_samples_list.txt | cut -f2 | sort -u)
echo "sample_name: ${SAMPLE_ID}"

# Define the VCFs and TSVs directory
TARGETS_DIR=${2}

# Define the VCF file
VCF=$(ls ${TARGETS_DIR} | grep "${SAMPLE_ID}" | grep "vcf.gz")
echo "TARGETS_DIR/VCF: ${TARGETS_DIR}/${VCF}"

# Define the TSV file
TSV=$(ls ${TARGETS_DIR} | grep "${SAMPLE_ID}" | grep "tsv.gz" | grep -v ".tbi")
echo "TARGETS_DIR/TSV: ${TARGETS_DIR}/${TSV}"

# Define the output directory
OUT_DIR=$(echo "${TARGETS_DIR}" | sed 's/ref_panels/genotype_likelihoods/g')
mkdir -p ${OUT_DIR}
echo "Output directory: ${OUT_DIR}"

# Define the output file
OUT=$(echo "${OUT_DIR}/${SAMPLE_ID}_${COV}_mLynPar1.2_ref_GL.vcf.gz")
echo "Output file: ${OUT}"

# Compute Genotypes Likelihoods
bcftools mpileup -f ${REFGEN} -I -E -a 'FORMAT/DP' -T  ${TARGETS_DIR}/${VCF} ${BAM} -Ou | bcftools call -Aim -C alleles -T ${TARGETS_DIR}/${TSV} -Oz -o ${OUT}
bcftools index -f ${OUT}


