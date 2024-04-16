#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/phasing/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/phasing/slurm-%j.out
#SBATCH -t 1-00:00:00
#SBATCH -c 5
#SBATCH --mem=20GB

# This script will run WhatsHap for phase set tagging on the selected chromosome.

# Usage of the script: chr_vcf_whatshap.sh <input_vcf> <bams_directory>

# charge the WhatsHap module:
module load cesga/2020 whatshap/1.1

# define input vcf
input_vcf=${1}

# define output directory
vcf_dir=$(dirname ${input_vcf})

# define input vcf basename
vcf_basename=$(basename ${input_vcf} .vcf)

# define reference genome
ref=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.scaffolds.revcomp.scaffolds.fa

# define bams directory
bams_dir=${2}

# WhatsHap for PS
whatshap phase \
 -o ${vcf_dir}/${vcf_basename}_ps.vcf \
 --tag=PS \
 --reference=${ref} \
 ${input_vcf} \
 $(ls ${bams_dir}/*_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.bam)
