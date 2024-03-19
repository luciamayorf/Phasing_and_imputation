# Phasing_and_imputation

In this repository, I will perform the phasing and imputation of data sequenced at low-coverage using a reference panel of 50 individuals sequenced at ~25X. I will start by phaing the VCF of the high-coverage sequenced invididuals.

## 1. Reference panel VCF phasing

First of all, we need to phase the VCF of the reference panel (high coverage sequenced individuals). For that, I will used the combined [WhatsHap v1.1](https://whatshap.readthedocs.io/en/latest/index.html) and [SHAPEIT4](https://odelaneau.github.io/shapeit4/) approach, as in [Enrico's Lynxtrongression repository](https://github.com/Enricobazzi/Lynxtrogression)

The following pipeline first uses WhatsHap to create phase sets from individual read and population data. The output of WhatsHap is then passed to SHAPEIT4, that will infer the haplotypes of each sample for each chromosome.

### Splitting the VCF into chromosomes

To divide my VCF into single chromosome VCFs I ran a custom bash script [chr_vcf_split.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/phasing/chr_vcf_split.sh) <input_vcf> <chr_bed>. The chromosomes I decided to keep are the larger ones: 18 autosomes and the X and Y chromosomes.

```bash
sbatch  -c 5 --mem=5GB -t 00:15:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/chr_vcf_split.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.vcf /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.big_chromosomes.bed
```

### Generate genetic map

To run SHAPEIT4 I also need to provide a genetic map for the SNPs to phase. As we don't have one, we will manually generate a genetic map by multiplying the physical distance in bp between SNPs and genome wide average recombination rate, which is 1.9 cM/Mbp. By cumulatively summing the multiplication of the physical distance from previous the SNP by 0.0000019, we obtain the cM value of each SNP. This approximation is not ideal but it's the only way we can provide a map. To calculate this I wrote a custom script [make_chr_gmap.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/phasing/make_chr_gmap.sh) <input_vcf> <chr_bed> which will output a gmap table for each chromosome, made of 3 columns: position, chromosome, cM (format useful for SHAPEIT4).

```bash
sbatch  -c 5 --mem=5GB -t 00:30:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/make_chr_gmap.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.vcf /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.big_chromosomes.bed
```

### Generate Phase sets with WhatsHap

For more precise phasing, we first run the software WhatsHap using the --tag=PS (see link).

Phase sets were generated from the VCF of each chromosome of each population by running in parallel a custom script pop_chr_vcf_whatshap.sh.


