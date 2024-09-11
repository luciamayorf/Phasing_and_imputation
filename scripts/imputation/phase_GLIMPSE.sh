#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/slurm-%j.out

# This script performs the multi-target phasing and imputation with GLIMPSE (Second step).
# The idea is to launch one job per chromosome. Usage of the script: 
# for vcf_gl in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_all/genotypes_likelihoods/epil_medcov_GL_merged.*.vcf.gz); do
#   sbatch -t 01:00:00 --mem 3GB phase_GLIMPSE.sh ${vcf_gl} <input_VCFref>
# done

# Load the GLIMPSE1.1 module
module load cesga/2020 gcc/system glimpse/1.1.1

# Define the input VCF with Genotype Likelihoods
VCFgl=${1}
echo "Input VCF (with GLs, to be imputed): ${VCFgl}"

# Define the input reference panel VCF file
VCF_REF=${2}
echo "Reference VCF: ${VCF_REF}"

# Define the imputation prefix (batch of samples imputed ID)
ID_IMP=$(basename ${VCFgl} .vcf.gz | sed 's/_GL_merged.mLynPar1.2_.*//')
echo "Batch imputation prefix: ${ID_IMP}"

# Define the chromosome
CHR=$(basename ${VCFgl} .vcf.gz | sed "s/${ID_IMP}_GL_merged.mLynPar1.2_//")
echo "Chromosome to be phased: ${CHR}"

# Define the output directory
OUT_DIR=$(dirname ${VCFgl} | sed 's/genotypes_likelihoods/GLIMPSE_phase/g')/${CHR}
mkdir -p ${OUT_DIR}
echo "Output directory of phasing and imputation step: ${OUT_DIR}"

# Define the chunk file
CHUNK_FILE=$(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/GLIMPSE_chunk/*_chunks.txt | grep "${CHR}")
echo "Chunk file: ${CHUNK_FILE}"


# Loop through each chunk file and perform phasing and imputation
while IFS="" read -r LINE || [ -n "$LINE" ]; do
    printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)

    # Define the variables
    IRG=$(echo $LINE | cut -d" " -f3)
    ORG=$(echo $LINE | cut -d" " -f4)
    OUT=${OUT_DIR}/${ID_IMP}_${IRG}_${ORG}.bcf
    echo "Input region: ${IRG}"
    echo "Output region: ${ORG}"
    echo "Output file: ${OUT}"

    # Perform phasing and imputation and index the output
    GLIMPSE_phase --input ${VCFgl} --reference ${VCF_REF} \
    --map /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/gmaps/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_originalnames_${CHR}.gmap \
    --input-region ${IRG} --output-region ${ORG} --output ${OUT}
    bcftools index -f ${OUT}

done < ${CHUNK_FILE}
