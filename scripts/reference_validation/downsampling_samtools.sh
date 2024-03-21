#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/downsampling/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/downsampling/slurm-%j.out

# This script performs the BAM downsampling to the 3 coverages displayed in a table with the structure: 
# <sample> <sample_coverage> <downsampling_factor1> <downsampling_factor2> <downsampling_factor3>

# Usage of the script: 
# for input_bam in $(ls /path/to/bams/*_.bam); do 
# job_id=$(sbatch -c 5 --mem=10GB -t 00:30:00 ./downsampling_samtools.sh ${input_bam} <ds_table.txt> | awk '{print $4}')
#   echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/downsampling/job_ids_downsampling_samtools.txt
# done 

# Load the samtools module:
module load cesga/2020 samtools/1.19

# Define input directory:
input_bam=${1}

# Define bams directory:
bam_dir=$(dirname ${input_bam})

# Define the bams basename:
bam_basename=$(basename ${input_bam} .bam)

# Define the downsampling factors table:
ds_list=${2}

# Define the sample name coverages: CAREFUL WITH THE CUT COMMAND, it depends on the sample name format in the bam file
sample=$(echo "${bam_basename}" | cut -f1,2 -d'_')

# Define the downsampling factors:
target_2=$(grep "${sample}" ${ds_list} | cut -f5)q
target_1=$(grep "${sample}" ${ds_list} | cut -f4)
target_0_5=$(grep "${sample}" ${ds_list} | cut -f3)



#### Subsampling the bam files to the target coverages

## Subsample to 2x
echo "Subsampling ${bam_basename} to 2x"
samtools view -@ 20 -bh -s ${target_2} ${input_bam} > ${bam_dir}/downsampling/${bam_basename}.2x.bam
samtools index ${bam_dir}/downsampling/${bam_basename}.2x.bam

## Subsample to 1x
echo "Subsampling ${bam_basename} to 1x"                                                             
samtools view -@ 20 -bh -s ${target_1} ${input_bam} > ${bam_dir}/downsampling/${bam_basename}.1x.bam
samtools index ${bam_dir}/downsampling/${bam_basename}.1x.bam

## Subsample to 0.5x
echo "Subsampling ${bam_basename} to 0.5x"                                                             
samtools view -@ 20 -bh -s ${target_0_5} ${input_bam} > ${bam_dir}/downsampling/${bam_basename}.0_5x.bam
samtools index ${bam_dir}/downsampling/${bam_basename}.0_5x.bam

