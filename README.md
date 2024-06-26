# Phasing_and_imputation

In this repository, I will perform the phasing and imputation of data sequenced at low-coverage using a reference panel of 50 individuals sequenced at ~25X. 


## 1. Reference panel VCF phasing

First of all, we need to phase the VCF of the reference panel (50 high coverage sequenced individuals). For that, I will used the combined [WhatsHap v1.1](https://whatshap.readthedocs.io/en/latest/index.html) and [SHAPEIT4](https://odelaneau.github.io/shapeit4/) approach, as in [Enrico's Lynxtrongression repository](https://github.com/Enricobazzi/Lynxtrogression).

The following pipeline first uses WhatsHap to create phase sets from individual read and population data. The output of WhatsHap is then passed to SHAPEIT4, that will infer the haplotypes of each sample for each chromosome.

I need to change the name of the samples again in the VCF because they need to match the read groups from the BAM files for WhatsHap. I should've changed them before doing the alignment to avoid this extra steps (carefull, the order of the samples was alphabetically changed in the step 2 of the filtering):

```bash
module load samtools
bcftools reheader -s <(sort -k2 /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/FASTQ_files/novogene_lp_sept2023/fastq_samples_list.txt | cut -f1,2 -d'_' | uniq) -o c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_originalnames.vcf c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.vcf
```

### Splitting the VCF into chromosomes

To divide my VCF into single chromosome VCFs I ran a custom bash script [chr_vcf_split.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/phasing/chr_vcf_split.sh) <input_vcf> <chr_bed>. The chromosomes I decided to keep are the larger ones: 18 autosomes and the X chromosome.

```bash
sbatch --mem=5GB -t 00:15:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/chr_vcf_split.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_originalnames.vcf /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.big_chromosomes.bed
```

### Generate genetic map

To run SHAPEIT4 I also need to provide a genetic map for the SNPs to phase. As we don't have one, we will manually generate a genetic map by multiplying the physical distance in bp between SNPs and genome wide average recombination rate, which is 1.9 cM/Mbp. By cumulatively summing the multiplication of the physical distance from previous the SNP by 0.0000019, we obtain the cM value of each SNP. This approximation is not ideal but it's the only way we can provide a map. To calculate this I wrote a custom script [make_chr_gmap.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/phasing/make_chr_gmap.sh) <input_vcf> <chr_bed> which will output a gmap table for each chromosome, made of 3 columns: position, chromosome, cM (format useful for SHAPEIT4).

```bash
sbatch --mem=5GB -t 00:30:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/make_chr_gmap.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_originalnames.vcf /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.big_chromosomes.bed
```

### Generate Phase sets with WhatsHap

For a more precise phasing, we first run the software WhatsHap using the --tag=PS (see link). Phase sets were generated from the VCF of each chromosome of each population by running in parallel a custom script [chr_vcf_whatshap.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/phasing/chr_vcf_whatshap.sh) <input_vcf> <bams_directory>.

```bash
for input_vcf in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/chr_vcfs/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_originalnames_*.vcf); do 
  job_id=$(sbatch /home/csic/eye/lmf/scripts/Phasing_and_imputation/chr_vcf_whatshap.sh ${input_vcf} /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23 | awk '{print $4}')
    echo "${job_id} ${input_vcf}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/phasing/job_ids_chr_vcf_whatshap_novogene_lp_sept2023.txt
done
```

### Phase using SHAPEIT4

The data is now ready to be phased using SHAPEIT4. To do so in parallel, I used a custom made script [chr_vcf_shapeit.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/phasing/chr_vcf_shapeit.sh) <input_vcf>, that runs SHAPEIT4 for each chromosome, zipping the file and indexing it (necessary to run SHAPEIT4). MCMC iterations were set to "10b,1p,1b,1p,1b,1p,1b,1p,10m" as suggested by the SHAPEIT4 manual.

```bash
for input_vcf in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/chr_vcfs/*_ps.vcf); do 
  job_id=$(sbatch /home/csic/eye/lmf/scripts/Phasing_and_imputation/chr_vcf_shapeit.sh ${input_vcf} | awk '{print $4}')
    echo "${job_id} ${input_vcf}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/phasing/job_ids_chr_vcf_shapeit_novogene_lp_sept2023.txt
done 
```
We decide to keep the imputed genotypes. In the end, we are interested in providing a set of reference haplotypes. Therefore, missing genotypes only represent alternative haplotypes that would not be present in the reference panel, and would therefore be impossible to impute. In conclusion, the final set of reference haplotypes will be the similar with or without the missing genotypes. If at some point we are interested in setting those genotypes as missing for downstream analysis, check Enrico's [gt_masker_pop_chr_vcf.sh](https://github.com/Enricobazzi/Lynxtrogression/blob/main/scripts/phasing/gt_masker_pop_chr_vcf.sh) script.


### Phased VCF merging

Now we will merge all the chromosomes files to obtain a final VCF with all the phased variants (a total of 1324598 SNPs).
```bash
module load bcftools
vcf-concat $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/chr_vcfs/phasing/ | grep -v "ChrY") > ./../../c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_originalnames.phased.vcf
```

I change again the names of the VCF file so that they have the correct names and index it:
```bash
bcftools reheader -s <(sort -k2 /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/FASTQ_files/novogene_lp_sept2023/fastq_samples_list.txt | cut -f2 | uniq) -o c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.vcf c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_originalnames.phased.vcf

bgzip c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.vcf > c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.vcf.gz
tabix -p vcf c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.vcf
```

We need to retrieve the old header, where most information was stripped away during phasing, otherwise I can get some errors when running different tools. To "fix" the header we simply copy the pre-phased header, add the few new fields added during phasing, and finally add the phased part of the table.
```bash
# take pre-phased header
grep "##" c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.vcf > c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.fixed.vcf

# add info and format added in phasing
grep -E "##INFO|##FORMAT" c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.vcf >> c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.fixed.vcf

# add phased vcf table
grep -v "##" c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.vcf >> c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.fixed.vcf
```


---

## 2. Imputation with GLIMPSE

### Computation of genotypes likelihoods

Before proceeding with imputation, I need to calculate the genotype likelihoods (GLs) of the targeted variants in our samples. For that, I will use [BCFtools](https://samtools.github.io/bcftools/bcftools.html#mpileup) mpileup and call, following [GLIMPSE manual](https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html#run_preliminaries) recommendations.

For that, I first need to generate and index a TSV file of the reference panel VCF. 
```bash
module load samtools/1.19

# Generate the TSV file
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.vcf.gz | bgzip -c > c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.tsv.gz

# Index the TSV file
tabix -s1 -b2 -e2 c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.tsv.gz
```

Then I can compute the GLs using BCFtools with the custom script [gl_bcftools.sh]() <input_bam> <reference_vcf> <output_directory>.
```{bash}
for input_bam in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/pool_epil_all/*_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.bam); do
  job_id=$(sbatch --mem=2GB -t 00:15:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/gl_bcftools.sh ${input_bam} /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.vcf.gz /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_all/genotypes_likelihoods | awk '{print $4}')
  echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/job_ids_gl_bcftools.txt
done
```



