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

Then I can compute the GLs using BCFtools with the custom script [gl_bcftools.sh]() <input_bam> <reference_vcf> <output_directory>. The output of the script is one VCF per sample that contains the GLs in the PL field of all the targeted variants.

In this case we are going to impute three different sets of samples: the low-coverage epileptic set, the medium coverage samples and the genome project ones. I will compute their GLs separately.
```{bash}
# For the epileptic set
for input_bam in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/pool_epil_all/*_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.bam); do
  job_id=$(sbatch --mem=2GB -t 00:15:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/gl_bcftools.sh ${input_bam} /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.vcf.gz /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_all/genotypes_likelihoods | awk '{print $4}')
  echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/job_ids_gl_bcftools.txt
done

# For the medium coverage set
for input_bam in $(ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/old_sequences/*_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.bam | grep -f <(cut -f15 -d'/' /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_gvcfs/old_sequences/list_gvcfs_old_sequences_medcov.txt | cut -f1-4 -d'_')); do
  job_id=$(sbatch --mem=3GB -t 00:20:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/gl_bcftools.sh ${input_bam} /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.fixed.vcf.gz /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/old_sequences/genotypes_likelihoods/medcov | awk '{print $4}')
  echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/job_ids_gl_bcftools.txt
done

# For the genome project samples
for input_bam in $(ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/old_sequences/*_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.bam | grep -f <(cut -f15 -d'/' /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_gvcfs/old_sequences/list_gvcfs_old_sequences_highcov.txt | cut -f1-4 -d'_')); do
  job_id=$(sbatch --mem=2GB -t 01:00:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/gl_bcftools.sh ${input_bam} /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.fixed.vcf.gz /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/old_sequences/genotypes_likelihoods/genome_project | awk '{print $4}')
  echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/job_ids_gl_bcftools.txt
done
```

As GLIMPSE version 1.1 performs a multi-target imputation, the GLs of the different samples must be merged together to generate a single VCF. 

```bash
# generate the list of the vcfs to be merged:
ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_medcov/genotypes_likelihoods/*_GL.vcf.gz > /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_medcov_gp/list_epil_medcov_gp.txt
ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/old_sequences/genotypes_likelihoods/medcov/*_GL.vcf.gz >> /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_medcov_gp/list_epil_medcov_gp.txt
ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/old_sequences/genotypes_likelihoods/genome_project/*_GL.vcf.gz >> /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_medcov_gp/list_epil_medcov_gp.txt


# merge, separate by chromosomes and index the vcfs containing the GLs of the gp, the medium coverage and the epil samples:
module load samtools
for chr in $(cut -f1 /mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.big_chromosomes.bed); do
  bcftools merge -m none -r ${chr} -Oz -o /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_medcov_gp/genotypes_likelihoods/epil_medcov_gp_GL_merged.${chr}.vcf.gz -l /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_medcov_gp/list_epil_medcov_gp.txt
  bcftools index -f /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_medcov_gp/genotypes_likelihoods/epil_medcov_gp_GL_merged.${chr}.vcf.gz
done
```

### Chunks definition

This step only needs to be run once for the reference panel. I generated a custom script [chunks_GLIMPSE.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/imputation/chunks_GLIMPSE.sh). 
```bash
for chr in $(cut -f1 /mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.big_chromosomes.bed); do
  job_id=$(sbatch -t 00:10:00 --mem 2GB /home/csic/eye/lmf/scripts/Phasing_and_imputation/chunks_GLIMPSE.sh ${chr} /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.fixed.vcf.gz| awk '{print $4}')
  echo "${job_id} ${chr}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/job_ids_chunks.txt
done
```

### Phasing

I use the custom script [phase_GLIMPSE.sh](https://github.com/luciamayorf/Phasing_and_imputation/tree/main/scripts/imputation). This step imputates the genotypes, generating one bcf file per chunk.

```{bash}
for vcf_gl in $(ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_medcov_gp/genotypes_likelihoods/epil_medcov_gp_GL_merged.*.vcf.gz | grep -v "ChrX"); do
  job_id=$(sbatch -t 01:00:00 --mem 2GB /home/csic/eye/lmf/scripts/Phasing_and_imputation/phase_GLIMPSE.sh ${vcf_gl} /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.phased.fixed.vcf.gz | awk '{print $4}')
  echo "${job_id} ${chr}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/job_ids_phase.txt
done
```
### Ligation

Later, the all the bcf files generated need to be ligated for each chromosome, using the script [ligate_GLIMPSE.sh](https://github.com/luciamayorf/Phasing_and_imputation/tree/main/scripts/imputation)

```bash
for chr_dir in $(ls -d /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_medcov_gp/GLIMPSE_phase/*); do
  job_id=$(sbatch -t 00:20:00 --mem 2GB /home/csic/eye/lmf/scripts/Phasing_and_imputation/ligate_GLIMPSE.sh ${chr_dir} epil_medcov_gp | awk '{print $4}')
  echo "${job_id} ${chr_dir}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/job_ids_ligate.txt
done
```

### Sampling

This steps phases the genotypes, following the script [sample_GLIMPSE.sh](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/imputation/sample_GLIMPSE.sh)
```bash
for bcf in $(ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_medcov_gp/GLIMPSE_ligate/*.bcf); do
  job_id=$(sbatch -t 00:20:00 --mem 1GB /home/csic/eye/lmf/scripts/Phasing_and_imputation/sample_GLIMPSE.sh ${bcf} | awk '{print $4}')
  echo "${job_id} ${bcf}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/job_ids_sample.txt
done
```

### Concatenating autosomes files

As a final step, we want to merge all the individual autosomes files into a single file, for both the sampled and ligated BCFs. 

```bash
# concatenate and index the ligated files
bcftools concat $(find /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_medcov_gp/GLIMPSE_ligate -name "*.bcf" | grep -v "ChrX") -O b -o /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_medcov_gp/GLIMPSE_ligate/epil_medcov_gp_autosomes.bcf
bcftools index -f /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_medcov_gp/GLIMPSE_ligate/epil_medcov_gp_autosomes.bcf

# concatenate and index the sampled files
bcftools concat $(find /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_medcov_gp/GLIMPSE_sample -name "*.bcf" | grep -v "ChrX") -O b -o /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_medcov_gp/GLIMPSE_sample/epil_medcov_gp_autosomes_phased.bcf
bcftools index -f /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_medcov_gp/GLIMPSE_sample/epil_medcov_gp_autosomes_phased.bcf
```
---

## 3. Imputation QC

### Global QC

We will analyse the values of the INFO field and the maximum genotype probability (maxGP - each genotype has 3 probabilities: probability of the genotype being AA, AB or BB). 

```bash
# INFO file
bcftools query -f '%INFO/INFO\n' epil_medcov_gp_autosomes.bcf > imputation_QC/info_epil_medcov_gp.txt

# MaxGP file (number of lines is the n_SNPs*n_samples)
bcftools query -f '[%SAMPLE=%GP\n]' epil_medcov_gp_autosomes.bcf | awk -F'=' '{split($2,a,","); print $1, (a[1]>a[2])?(a[1]>a[3]?a[1]:a[3]):(a[2]>a[3]?a[2]:a[3])}' > imputation_QC/maxGP_epil_medcov_gp.txt
```

The script [global_INFO_maxGP_plots.R](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/imputation/global_INFO_maxGP_plots.R) plots these two distributions
```bash
Rscript /home/csic/eye/lmf/scripts/Phasing_and_imputation/imputation_QC/global_info_maxGP_plots.R /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_medcov_gp/GLIMPSE_ligate/imputation_QC epil_medcov_gp
```

### Individual QC

We are going to analyze the individual patterns of the maxGP, by individual samples and by sequencing batchs (pool_epil, gp and medcov), by first diving the BCFs:

```bash
# To separate the BCF per batch
bcftools view -S <(cut -f15 -d'/' /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_gvcfs/old_sequences/list_gvcfs_old_sequences_medcov.txt | cut -f1-4 -d'_') epil_medcov_autosomes.bcf > medcov_autosomes.vcf
bcftools view -S <(ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_medcov/genotypes_likelihoods/c_lp_*.vcf.gz | awk -F'/' '{print $NF}' | cut -f1-4 -d'_') epil_medcov_autosomes.bcf > epil_autosomes.vcf
bcftools view -S <(cut -f15 -d'/' /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_gvcfs/old_sequences/list_gvcfs_old_sequences_highcov.txt | cut -f1-4 -d'_') epil_medcov_gp_autosomes.bcf > gp_autosomes.vcf

 ## CAREFUL, AF and INFO tags are not recalculated, but I only want to keep the maxGPs.

# To obtain the maxGP tables
for vcf in $(ls *_autosomes.vcf); do
    BATCH=$(echo ${vcf} | sed 's/_autosomes.vcf//')
    echo "Processing ${vcf}"
    bcftools query -f '[%SAMPLE=%GP\n]' ${vcf} | awk -F'=' '{split($2,a,","); print $1, (a[1]>a[2])?(a[1]>a[3]?a[1]:a[3]):(a[2]>a[3]?a[2]:a[3])}' > imputation_QC/maxGP_per_sample_${BATCH}.txt
done
```
These tables contain one row PER GENOTYPE, containing the sample name and the max GP for that genotype. As a result, samples appear as many times as genotypes they have.

Plots of these distributions can be obtained with [maxGP_plot.R](https://github.com/luciamayorf/Phasing_and_imputation/blob/main/scripts/imputation/maxGP_plot.R)

```bash
for vcf in $(ls *_autosomes.vcf); do
    BATCH=$(echo ${vcf} | sed 's/_autosomes.vcf//')
    echo "Generating plots of ${BATCH}"
    Rscript /home/csic/eye/lmf/scripts/Phasing_and_imputation/imputation_QC/maxGP_plot.R /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/imputation_GLIMPSE/pool_epil_medcov_gp/GLIMPSE_ligate/imputation_QC maxGP_per_sample_${BATCH}
done
```

We also want to obtain a general metric per individual, based on the maxGP. We decided to study the number of SNPs with a maxGP under 0.95 per sample (a way to measure how many genotypes would not be reliable). To get those tables:

```bash
for vcf in $(ls *_autosomes.vcf); do
    BATCH=$(echo ${vcf} | sed 's/_autosomes.vcf//')
    echo "Generating maxGP< 0.95 table of ${BATCH}"
    
    for sample in $(bcftools query -l ${vcf}); do
      MAXGP=$(grep "${sample}" imputation_QC/maxGP_per_sample_${BATCH}.txt | awk '$2 < 0.95' | wc -l)
      echo "${sample} ${MAXGP}" >> imputation_QC/maxGP095_per_sample_${BATCH}.txt
    done
    
done
```


