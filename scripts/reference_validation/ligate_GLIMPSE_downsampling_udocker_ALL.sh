#!/bin/bash

# To use other samples I would need to change the INPUT_DIR and the OUT_DIR

# Loop trough each directory containing the chunk bcf files
for sample_cov in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/ref_panel_validation/GLIMPSE_phase) ; do

    # Define the input directory
    INPUT_DIR=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/ref_panel_validation/GLIMPSE_phase/${sample_cov}
    echo "Input directory of ligating step: ${INPUT_DIR}"

    # Define the output directory
    OUT_DIR=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/ref_panel_validation/GLIMPSE_ligate/${sample_cov}
    mkdir -p ${OUT_DIR}
    echo "Output directory of ligating step: ${OUT_DIR}"
    
    # Loop through the chromosomes files and ligate them 
    for chr in $(cut -f1 /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.big_chromosomes.bed); do
        
        # Define list of BCF files to be ligated
        find ${INPUT_DIR} -name "*.bcf" | grep "${chr}" > ${INPUT_DIR}/${sample_cov}_${chr}_bcf_list.txt

        # Ligate the BCF files and index the output
        GLIMPSE_ligate_v1.1.1 --input ${INPUT_DIR}/${sample_cov}_${chr}_bcf_list.txt --output ${OUT_DIR}/${sample_cov}_${chr}_merged.bcf 
        bcftools index -f ${OUT_DIR}/${sample_cov}_${chr}_merged.bcf 

        # Delete the list of BCF files
        rm ${INPUT_DIR}/${sample_cov}_${chr}_bcf_list.txt
        
    done
done
