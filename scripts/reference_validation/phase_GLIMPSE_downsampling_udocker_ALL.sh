#!/bin/bash

# To use other samples I would need to change the INPUT_DIR and the OUT_DIR

# Define the sample input directory containing all the chunk txt files
INPUT_DIR=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/ref_panel_validation/GLIMPSE_chunk
echo "Input directory is: ${INPUT_DIR}"

for input_vcf in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/ref_panel_validation/genotype_likelihoods/*_ref_GL.vcf.gz) ; do
    
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
    VCF_REF=$(ls ${VCF_REF_DIR} | grep "${SAMPLE}" | grep ".vcf.gz" | grep -v ".csi")
    echo "Reference VCF: ${VCF_REF_DIR}/${VCF_REF}"

    # Define the first step output directory
    OUT_DIR=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/ref_panel_validation/GLIMPSE_phase/${SAMPLE_COV}
    mkdir -p ${OUT_DIR}
    echo "Output directory of phasing step: ${OUT_DIR}"

    ### 2. PHASING AND IMPUTATION ###

    # Loop to each chunk file and perform phasing and imputation
    for chunk_file in $(ls ${INPUT_DIR}/${SAMPLE_COV}); do
        while IFS="" read -r LINE || [ -n "$LINE" ]; do
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)

            # Define the variables
            IRG=$(echo $LINE | cut -d" " -f3)
            ORG=$(echo $LINE | cut -d" " -f4)
            OUT=${OUT_DIR}/${SAMPLE_COV}_${IRG}_${ORG}.bcf
            CHR=$(echo $LINE | cut -d" " -f2 | sed 's/mLynPar1.2_//g')

            # Perform phasing and imputation and index the output
            GLIMPSE_phase_v1.1.1 --input ${VCFgl} --reference ${VCF_REF_DIR}/${VCF_REF} \
            --map /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/gmaps/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_originalnames_${CHR}.gmap \
            --input-region ${IRG} --output-region ${ORG} --output ${OUT}
            bcftools index -f ${OUT}
        done < ${INPUT_DIR}/${SAMPLE_COV}/${chunk_file}
    done
done