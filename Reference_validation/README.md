# Leave-one-out downsampling and imputation

In this section we will perform the validation of the reference to estimate imputations accuracy. We will follow the leave-one-out approach to do so.


## BAMs downsampling

First of all, we need to downsample the BAM files of the high coverage individuals to replicate the coverage we have for the low-coverage data. We will target 3 coverages: 0.5X, 1X and 2X.

### Calculating the downsampling factor

To do so, we first need to calculate the proportion of reads that need to be subsampled for each sample to get the target coverages. For that, I will first calculate the coverage per sample with [mosdepth](https://github.com/brentp/mosdepth) using the script [mean_cov_mosdepth.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/reference_validation/mean_cov_mosdepth.sh) <input_bam> <output_directory>:

```bash
for input_bam in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/*_sorted_rg_merged_sorted_rmdup_indelrealigner.bam); do 
  job_id=$(sbatch -c 16 --mem=10GB -t 00:15:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/ref_panel_validation/mean_cov_mosdepth.sh ${input_bam} /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/mosdepth | awk '{print $4}')
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
  job_id=$(sbatch -c 5 --mem=10GB -t 00:45:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/ref_panel_validation/downsampling_samtools.sh ${input_bam} /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/mosdepth/downsampling_factor_table.txt | awk '{print $4}')
    echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/downsampling/job_ids_downsampling_samtools.txt
done
```

#### QC of downsampled BAM

Now I need to check that the coverage of the downsampled files is similar to our target. I will also perform a qualimap bamqc to see the distribution of reads along the genome and compare if it similar to our lcWGS data.

I run mosdepth to get the coverage with the same script as before ([mean_cov_mosdepth.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/reference_validation/mean_cov_mosdepth.sh)):

```bash
for input_bam in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/downsampling/*.bam); do 
  job_id=$(sbatch -c 5 --mem8GB -t 00:10:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/ref_panel_validation/mean_cov_mosdepth.sh ${input_bam} /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/downsampling/mosdepth | awk '{print $4}')
    echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/mosdepth/job_ids_mean_cov_mosdepth.txt
done
```

I then perform a bam quality control with qualimap by running the script [bams_qualimap.sh](https://github.com/luciamayorf/Data_preprocessing_alignment_v2/blob/main/scripts/bams_qualimap.sh)
```bash
for input_bam in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/downsampling/*.bam); do
  job_id=$(sbatch -c 5 --mem=10GB -t 00:20:00 /home/csic/eye/lmf/scripts/Data_preprocessing_alignment/bams_qualimap.sh ${input_bam} | awk '{print $4}')
  echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/qualimap/job_ids_qualimap_dowmsampling.txt
done
```

Finally, I run multiqc to visualize the results.

---

## Computation of genotype likelihoods

GLIMPSE requires input data to take the form of Genotype Likelihoods (GLs). GLs need to be computed at all target individuals and all variant sites present in the reference panel of haplotypes used for the imputation. We used [BCFtools](https://samtools.github.io/bcftools/bcftools.html#mpileup) mpileup and call, as it was done in the [GLIMPSE manual](https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html#run_preliminaries). 

For that, I first need to generate different VCF files, one for each sample extracted (that will serve as a reference panel for its downsampled BAM). I will run the custom script [sample_removal_vcf_downsampling.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/reference_validation/sample_removal_vcf_downsampling.sh) <input_vcf>, which also generates the TSV files necessary for bcftools to calculate the genotypes likelihoods.
```bash
sbatch --mem=1GB -t 00:20:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/ref_panel_validation/sample_removal_vcf_downsampling.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.vcf.gz

# I moved all the generated files to a new directory
mv * ref_panels/
```

Now we can compute the GLs using BCFtools with the custom script [gl_bcftools_downsampled.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/reference_validation/gl_bcftools_downsampled.sh) <input_bam> <vcf_tsv_directory>. I run each coverage separately to avoid exceeding CESGA jobs launching limits.

```bash
for input_bam in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/downsampling/*0_5x.bam); do
  job_id=$(sbatch --mem=2GB -t 00:10:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/ref_panel_validation/gl_bcftools_downsampled.sh ${input_bam} /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/ref_panel_validation/ref_panels | awk '{print $4}')
  echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/downsampling/job_ids_gl_bcftools_downsampled.txt
done

for input_bam in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/downsampling/*1x.bam); do
  job_id=$(sbatch --mem=2GB -t 00:10:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/ref_panel_validation/gl_bcftools_downsampled.sh ${input_bam} /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/ref_panel_validation/ref_panels | awk '{print $4}')
  echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/downsampling/job_ids_gl_bcftools_downsampled.txt
done

for input_bam in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/downsampling/*2x.bam); do
  job_id=$(sbatch --mem=2GB -t 00:15:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/ref_panel_validation/gl_bcftools_downsampled.sh ${input_bam} /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/ref_panel_validation/ref_panels | awk '{print $4}')
  echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/downsampling/job_ids_gl_bcftools_downsampled.txt
done
```

---

## Reference panel imputation

Once the BAMs of the samples in the reference panel are downsampled to different target coverages and their GLs are computed, we can start using GLIMPSE to impute them.

Note: all these scripts were executed when GLIMPSEv1.1 was still not installed in CESGA ft3. I had to download the docker image and run the scripts a little bit differently. I need to first opne the Docker image, using the script [launch_scriptGLIMPSE_downsampling.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/reference_validation/launch_scriptGLIMPSE_downsampling.sh), followed by the script I want to run afterwards.

Following [GLIMPSE1 tutorial](https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html#run_preliminaries), I run one script per stage of the whole imputation process.
 

### Chunks definition

I tried generating one script per sample, but the image was "broken" when running simultaneouly. In the end, I decide to generate a script with a loop that does them all at the same time (not ideal, but it worked and it didn't take that long). To define the chunks, I run the custom script [chunks_GLIMPSE_downsampling_udocker_ALL.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/reference_validation/launch_scriptGLIMPSE_downsampling.sh) inside the Docker image.

 ```{bash}
sbatch -t 00:30:00 --mem 3GB /home/csic/eye/lmf/scripts/Phasing_and_imputation/ref_panel_validation/launch_scriptGLIMPSE_downsampling.sh \
/home/csic/eye/lmf/scripts/Phasing_and_imputation/ref_panel_validation/chunks_GLIMPSE_downsampling_udocker_ALL.sh   
```
This steps generates one TXT file per chunk.

 ### Phasing

 Once the chunks are generated, next step consistf of performing the phasing and imputation of missing genotypes, running the script [phase_GLIMPSE_downsampling_udocker_ALL.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/reference_validation/phase_GLIMPSE_downsampling_udocker_ALL.sh).

 ```bash
sbatch -t 1-00:00:00 --mem 3GB 
/home/csic/eye/lmf/scripts/Phasing_and_imputation/ref_panel_validation/launch_scriptGLIMPSE_downsampling.sh \
/home/csic/eye/lmf/scripts/Phasing_and_imputation/ref_panel_validation/phase_GLIMPSE_downsampling_udocker_ALL.sh
````

### Ligation

After phasing and imputation, one BCF is generated per chunk. The next step consists of merging all the chunk files of the same chromosome with the custom script [ligate_GLIMPSE_downsampling_udocker_ALL.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/reference_validation/ligate_GLIMPSE_downsampling_udocker_ALL.sh).

 ```bash
sbatch -t 02:00:00 --mem 2GB
/home/csic/eye/lmf/scripts/Phasing_and_imputation/ref_panel_validation/launch_scriptGLIMPSE_downsampling.sh \
/home/csic/eye/lmf/scripts/Phasing_and_imputation/ref_panel_validation/ligate_GLIMPSE_downsampling_udocker_ALL.sh
````

### Sampling haplotypes

This step outputs phased genotypes. We run the script [sample_GLIMPSE_downsampling_udocker_ALL.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/reference_validation/sample_GLIMPSE_downsampling_udocker_ALL.sh).

 ```bash
sbatch -t 01:00:00 --mem 2GB
/home/csic/eye/lmf/scripts/Phasing_and_imputation/ref_panel_validation/launch_scriptGLIMPSE_downsampling.sh \
/home/csic/eye/lmf/scripts/Phasing_and_imputation/ref_panel_validation/sample_GLIMPSE_downsampling_udocker_ALL.sh
````

### Concordance
