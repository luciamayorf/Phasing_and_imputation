# Leave-one-out downsampling and imputation

In this section we will perform the validation of the reference to estimate imputations accuracy. We will follow the leave-one-out approach to do so.

## BAMs downsampling

First of all, we need to downsample the BAM files of the high coverage individuals to replicate the coverage we have for the low-coverage data. We will target 3 coverages: 0.5X, 1X and 2X.

To do so, we first need to calculate the proportion of reads that need to be subsampled for each sample to get the target coverages. For that, I will first calculate the coverage per sample with [mosdepth](https://github.com/brentp/mosdepth) using the script [mean_cov_mosdepth.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/reference_validation/mean_cov_mosdepth.sh)

```bash
for input_bam in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/*_sorted_rg_merged_sorted_rmdup_indelrealigner.bam); do 
  job_id=$(sbatch -c 16 --mem=10GB -t 00:15:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/mean_cov_mosdepth.sh ${input_bam} /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/mosdepth | awk '{print $4}')
    echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/mosdepth/job_ids_mean_cov_mosdepth.txt
done
```
