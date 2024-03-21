# Leave-one-out downsampling and imputation

In this section we will perform the validation of the reference to estimate imputations accuracy. We will follow the leave-one-out approach to do so.

## BAMs downsampling

First of all, we need to downsample the BAM files of the high coverage individuals to replicate the coverage we have for the low-coverage data. We will target 3 coverages: 0.5X, 1X and 2X.

To do so, we first need to calculate the proportion of reads that need to be subsampled for each sample to get the target coverages. For that, I will first calculate the coverage per sample with [mosdepth](https://github.com/brentp/mosdepth)
