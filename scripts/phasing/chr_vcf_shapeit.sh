#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/phasing/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/phasing/slurm-%j.out
#SBATCH -t 00:20:00
#SBATCH -c 4
#SBATCH --mem=10GB

# This script will run SHAPEIT4 for phasing the phase set VCFs of a specific chromosome 

# Usage of the script: chr_vcf_shapeit.sh <input_vcf>

# charge the SHAPEIT4 module:
module load cesga/2020 gcccore/system shapeit4/4.2.1
module load cesga/2020 bcftools/1.19

# define input vcf
input_vcf=${1}

# define input vcf directory
vcf_dir=$(dirname ${input_vcf})

# define and create output directory if it doesn't exits
output_dir=$(echo "${vcf_dir}/phasing") 
mkdir -p ${output_dir}

# define input vcf basename
vcf_basename=$(basename ${input_vcf} .vcf)

# define genetic map
gmap=$(echo "${vcf_dir}" | sed 's/chr_vcfs/gmaps/')/$(basename ${input_vcf} _ps.vcf).gmap

# define region to be phased (the whole chromosome)
chr=$(grep -v "#" ${input_vcf} | head -n 1 | cut -f1)

# The vcf needs to be zipped and indexed:  
echo "zipping ${input_vcf}"
bgzip ${input_vcf}
echo "indexing ${input_vcf}.gz"
bcftools index ${input_vcf}.gz 

# SHAPEIT4 command
shapeit4.2 \
 --input ${input_vcf}.gz \
 --map ${gmap} \
 --region ${chr} \
 --use-PS 0.0001 \
 --output ${output_dir}/${vcf_basename}_phased.vcf \
 --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m \
 --thread 4