#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/phasing/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/phasing/slurm-%j.out

# This script generates a generic genetic map using a generic recombination rate of 1.9cM/Mb.

# Usage of the script: sbatch -c 5 --mem=10GB -t 01:00:00 make_chr_gmap.sh <input_vcf> <chr_bed>

# define input vcf
input_vcf=${1}

# define output directory
vcf_dir=$(dirname ${input_vcf})
mkdir -p ${vcf_dir}/gmaps

# define input vcf basename
vcf_basename=$(basename ${input_vcf} .vcf)

# define list of chromosomes
chr_list=${2}

# generate genetic map
for chr in $(cut -f1 ${chr_list}) ; do
    echo "calculating genetic map of ${chr} from ${input_vcf}"
    chr_name=$(echo "${chr}" | cut -d'_' -f2-)
    grep -v "#" ${input_vcf} | grep -w ${chr} | cut -f1,2 |
    awk '{ print $2, $1 }' |
    awk {'if ( NR==1 ) print $1, $2, 0; else print $1, $2, $1-p, ($1-p)*0.0000019; p=$1'} |
    awk 'BEGIN{print "pos", "chr", "cM"} {sum+=$4} {print $1, $2, sum}' |
    tr ' ' '\t' > ${vcf_dir}/gmaps/${vcf_basename}_${chr_name}.gmap
done
 