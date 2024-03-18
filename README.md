# Phasing_and_imputation

In this repository, I will perform the phasing and imputation of data sequenced at low-coverage using a reference panel of 50 individuals sequenced at ~25X. I will start by phaing the VCF of the high-coverage sequenced invididuals.

## 1. Reference panel VCF phasing

First of all, we need to phase the VCF of the reference panel (high coverage sequenced individuals). For that, I will used the combined [WhatsHap v1.1](https://whatshap.readthedocs.io/en/latest/index.html)v1.1 and [SHAPEIT4](https://odelaneau.github.io/shapeit4/) approach, as in [Enrico's Lynxtrongression repository](https://github.com/Enricobazzi/Lynxtrogression)

The following pipeline first uses WhatsHap to create phase sets from individual read and population data. The output of WhatsHap is then passed to SHAPEIT4, that will infer the haplotypes of each sample for each chromosome.

### Splitting the VCF into chromosomes

To divide my VCF into single chromosome VCFs I ran a custom bash script [pop_chr_vcf_split.sh](). The chromosomes I decided to keep are the larger ones: 18 autosomes and the X and Y chromosomes.

```bash
sbatch  -c 10 --mem=20GB -t 01:00:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/chr_vcf_split.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20.vcf /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.big_chromosomes.bed
```


