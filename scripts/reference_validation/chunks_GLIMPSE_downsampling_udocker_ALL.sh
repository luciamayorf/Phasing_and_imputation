#!/bin/bash

# To use other samples I would need to change the input VCFs directory and the OUT_DIR

# Define the output directory
OUT_DIR=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/ref_panel_validation/GLIMPSE_chunk
mkdir -p ${OUT_DIR}
echo "Output directory: ${OUT_DIR}"

# Loop though each VCF file: 
for input_vcf in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/ref_panel_validation/genotype_likelihoods/*_mLynPar1.2_ref_GL.vcf.gz) ; do

    # Define input VCF to be imputed (with Genotype Likelihoods)
    VCFgl=${input_vcf}
    echo "Input VCF (with GLs, to be imputed): ${VCFgl}"

    # Define sample name
    SAMPLE_COV=$(basename ${VCFgl} _mLynPar1.2_ref_GL.vcf.gz)
    echo "Sample name with coverage: ${SAMPLE_COV}"
    SAMPLE=$(echo ${SAMPLE_COV} | cut -d'_' -f1-4)
    echo "Sample name: ${SAMPLE}"

    # Define reference panel directory
    VCF_REF_DIR=$(dirname ${VCFgl} | sed 's/genotype_likelihoods/ref_panels/g')
    echo "Reference panel directory: ${VCF_REF_DIR}"

    # Define the reference panel VCF 
    VCF_REF=$(ls ${VCF_REF_DIR} | grep "${SAMPLE}" | grep ".vcf.gz")
    echo "Reference VCF: ${VCF_REF_DIR}/${VCF_REF}"

    # Define a new directory inside the ouput directory to store all the chunks from the same sample
    OUT_DIR_SAMP=${OUT_DIR}/${SAMPLE_COV}
    mkdir -p ${OUT_DIR_SAMP}
    echo "Output directory for chunks: ${OUT_DIR_SAMP}"


    ### 1. DEFINE CHUNKS ###

    # Loop through each chromosome and split into chunks
    for chr in $(cut -f1 /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.big_chromosomes.bed); do
        GLIMPSE_chunk_v1.1.1 --input ${VCF_REF_DIR}/${VCF_REF} --region ${chr} --window-size 2000000 --buffer-size 200000 --output ${OUT_DIR_SAMP}/${SAMPLE_COV}_${chr}_chunks.txt
    done

done
