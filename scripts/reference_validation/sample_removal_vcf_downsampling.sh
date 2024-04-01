#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/downsampling/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/downsampling/slurm-%j.out

# For each sample of a vcf phased with SHAPEIT4, this script generates a new vcf extracting that sample extracts. 
# It then adapts its format to the one required by bcftools call (TSV) call by using tabix.

# Usage of the script: sbatch sample_removal_vcf_downsampling.sh <input_vcf>

# Load the samtools module:
module load cesga/2020 samtools/1.19

# Define input vcf
input_vcf=${1}

# Define input basename
input_basename=$(basename ${input_vcf} .vcf.gz | sed 's/c_lp_all_novogene_sept23_//')

# Define output vcf directory
output_vcf_dir=$(dirname ${input_vcf})/ref_panel_validation
mkdir -p ${output_vcf_dir}

# Define the samples to be removed
samples=$(bcftools query -l ${input_vcf})

# Loop over the samples and generate a new vcf extracting each sample
for sample in ${samples}; do

    # Extract the sample and generate a new VCF file
    bcftools view -s ^${sample} -Oz -o ${output_vcf_dir}/no_${sample}_${input_basename}.vcf.gz ${input_vcf}

    # Generate the TSV file  
    bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' ${output_vcf_dir}/no_${sample}_${input_basename}.vcf.gz \
    | bgzip -c > ${output_vcf_dir}/no_${sample}_${input_basename}.tsv.gz

    # Index the TSV file
    tabix -s1 -b2 -e2 ${output_vcf_dir}/no_${sample}_${input_basename}.tsv.gz
done









