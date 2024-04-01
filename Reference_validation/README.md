# Leave-one-out downsampling and imputation

In this section we will perform the validation of the reference to estimate imputations accuracy. We will follow the leave-one-out approach to do so.


## BAMs downsampling

First of all, we need to downsample the BAM files of the high coverage individuals to replicate the coverage we have for the low-coverage data. We will target 3 coverages: 0.5X, 1X and 2X.

### Calculating the downsampling factor

To do so, we first need to calculate the proportion of reads that need to be subsampled for each sample to get the target coverages. For that, I will first calculate the coverage per sample with [mosdepth](https://github.com/brentp/mosdepth) using the script [mean_cov_mosdepth.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/reference_validation/mean_cov_mosdepth.sh) <input_bam> <output_directory>:

```bash
for input_bam in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/*_sorted_rg_merged_sorted_rmdup_indelrealigner.bam); do 
  job_id=$(sbatch -c 16 --mem=10GB -t 00:15:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/mean_cov_mosdepth.sh ${input_bam} /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/mosdepth | awk '{print $4}')
    echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/mosdepth/job_ids_mean_cov_mosdepth.txt
done
```

Now I will generate a table with the mean coverage per sample and the downsampling factor for each of the target coverages.

```bash
for input_mosdepth in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/mosdepth/*.summary.txt); do
  cov=$(grep -w "total" ${input_mosdepth} | cut -f4);
  sample=$(basename ${input_mosdepth} _mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner_mosdepth.mosdepth.summary.txt)
  target_0_5=$(LC_ALL=C printf "%.4f" $(echo "scale=4; 0.5/${cov}" | bc -l))
  target_1=$(LC_ALL=C printf "%.4f" $(echo "scale=4; 1/${cov}" | bc -l))
  target_2=$(LC_ALL=C printf "%.4f" $(echo "scale=4; 2/${cov}" | bc -l))
  echo -e "${sample}\t${cov}\t${target_0_5}\t${target_1}\t${target_2}" >> /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/mosdepth/downsampling_factor_table.txt
done
```
This table presents a column with the sample name, its mean coverage, and the factors to subsample the BAM file to obtain our 3 target coverages.


### Downsampling the BAM files

We will generate a BAM file with the 3 different coverages (0.5X, 1X and 2X) for each sample using [samtools v1.19](https://www.htslib.org/doc/samtools-view.html). I will use the script [downsampling_samtools.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/reference_validation/downsampling_samtools.sh) <input_bam> <table_ds_factors>

```bash
for input_bam in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/*_sorted_rg_merged_sorted_rmdup_indelrealigner.bam); do 
  job_id=$(sbatch -c 5 --mem=10GB -t 00:45:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/downsampling_samtools.sh ${input_bam} /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/mosdepth/downsampling_factor_table.txt | awk '{print $4}')
    echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/downsampling/job_ids_downsampling_samtools.txt
done
```

#### QC of downsampled BAM

Now I need to check that the coverage of the downsampled files is similar to our target. I will also perform a qualimap bamqc to see the distribution of reads along the genome and compare if it similar to our lcWGS data.

I run mosdepth to get the coverage with the same script as before ([mean_cov_mosdepth.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/reference_validation/mean_cov_mosdepth.sh)):

```bash
for input_bam in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/downsampling/*.bam); do 
  job_id=$(sbatch -c 5 --mem8GB -t 00:10:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/mean_cov_mosdepth.sh ${input_bam} /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/downsampling/mosdepth | awk '{print $4}')
    echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/mosdepth/job_ids_mean_cov_mosdepth.txt
done
```

I then perform a bam quality control with qualimap by running the script [bams_qualimap.sh](https://github.com/luciamayorf/Data_preprocessing_alignment_v2/blob/main/scripts/bams_qualimap.sh)
```bash
for input_bam in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/downsampling/*.bam); do
  job_id=$(sbatch -c 5 --mem=10GB -t 00:20:00 /home/csic/eye/lmf/scripts/Data_preprocessing_alignment/bams_qualimap.sh ${input_bam} | awk '{print $4}')
  echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/qualimap/job_ids_qualimap_dpwmsampling.txt
done
```

Finally, I run multiqc to visualize the results.

---

## Computation of genotype likelihoods

GLIMPSE requires input data to take the form of Genotype Likelihoods (GLs). GLs need to be computed at all target individuals and all variant sites present in the reference panel of haplotypes used for the imputation. We used BCFtools, as it was done in the [GLIMPSE manual](https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html#run_preliminaries). 

For that, I first need to generate different VCF files, one for each sample extracted (that will serve as a reference panel for its downsampled BAM). I will run the custom script [sample_removal_vcf_downsampling.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/reference_validation/sample_removal_vcf_downsampling.sh) <input_vcf>, which also generates the TSV files necessary for bcftools to calculate the genotypes likelihoods.
```bash
sbatch -c 1 --mem=1GB -t 00:20:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/sample_removal_vcf_downsampling.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.vcf.gz
```
