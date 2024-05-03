#!/bin/bash

# To use other samples I would need to change the INPUT_DIR and the OUT_DIR

# Loop trough each directory containing the chunk bcf files
for sample_cov in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/ref_panel_validation/GLIMPSE_ligate) ; do

    # Define the input directory
    INPUT_DIR=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/ref_panel_validation/GLIMPSE_ligate/${sample_cov}
    echo "Input directory of sampling step: ${INPUT_DIR}"

    # Define the output directory
    OUT_DIR=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/ref_panel_validation/GLIMPSE_sample/${sample_cov}
    mkdir -p ${OUT_DIR}
    echo "Output directory of ligating step: ${OUT_DIR}"
    
    # Loop through the chromosomes files and ligate them 
    for input_bcf in $(ls ${INPUT_DIR}/*.bcf); do

        # Define output VCF name
        OUT_VCF=$(basename ${input_bcf} | sed 's/_merged.bcf/_phased.bcf/g')
        
        # Sample haplotypes from the chr BCF files and index the output
        GLIMPSE_sample_v1.1.1 --input ${input_bcf} --solve --output ${OUT_DIR}/${OUT_VCF} 
        bcftools index -f ${OUT_DIR}/${OUT_VCF}  
        
    done

    # concatenate all the final BCF files from the autosomes (excluding the X chromosome)
    bcftools concat $(find ${OUT_DIR}  -name "*.bcf" | grep -v "ChrX") -O b -o ${OUT_DIR}/${sample_cov}_phased_autosomes.bcf
    bcftools index -f ${OUT_DIR}/${sample_cov}_phased_autosomes.bcf
    
done

