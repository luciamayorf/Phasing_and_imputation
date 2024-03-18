#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/phasing/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/phasing/slurm-%j.out

# This script splits a VCF file into single chromosomes using bedtools intersect.

# Usage of the script: sbatch -c 5 --mem=10GB -t 01:00:00 chr_vcf_split.sh <input_vcf> <chr_list>

# charge the bedtools module:
module load cesga/2020 bedtools/2.31.0

# define input vcf
input_vcf=${1}

# define output directory
vcf_dir=$(dirname ${input_vcf})
mkdir -p ${vcf_dir}/chr_vcfs

# define input vcf basename
vcf_basename=$(basename ${input_vcf})

# define list of chromosomes
chr_list=${2}

# split VCF into single chrs
for chr in $(cut -f1 ${chr_list}[@]) ; do
    echo "extracting ${chr} from VCF"
    bedtools intersect -header -a ${input_vcf} -b <(grep -w ${chr} ${chr_list}) > ${vcf_dir}/chr_vcfs/${vcf_basename}_${chr}.vcf
done
