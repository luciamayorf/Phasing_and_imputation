# Leave-one-out downsampling and imputation

In this section we will perform the validation of the reference to estimate imputations accuracy. We will follow the leave-one-out approach to do so.


## BAMs downsampling

First of all, we need to downsample the BAM files of the high coverage individuals to replicate the coverage we have for the low-coverage data. We will target 3 coverages: 0.5X, 1X and 2X.

### Calculation of downsampling factor

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

### Downsampling the BAM files

We will generate a BAM file with the 3 different coverages (0.5X, 1X and 2X) for each sample using [samtools v1.19](https://www.htslib.org/doc/samtools-view.html).

