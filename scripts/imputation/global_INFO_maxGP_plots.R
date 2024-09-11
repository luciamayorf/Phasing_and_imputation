# This R script generates a pdf plot that representes the INFO field of a TXT file generated from a VCF with the command: 
# bcftools query -f '%INFO/INFO\n' <bcf> > batch_name

# It also generates a pdf plot that representes the maximum GP distribution of a TXT file generated from a VCF with the command: 
# bcftools query -f '[%GP\n]' medcov_autosomes.vcf | awk -F',' '{print ($1>$2)?($1>$3?$1:$3):($2>$3?$2:$3)}'

# Usage: Rscript global_info_maxGP_plots.R <input_path> <batch_name>
library(ggplot2)

# Get the path and base name from the command-line arguments
path <- commandArgs(trailingOnly = TRUE)[1]
base_name <- commandArgs(trailingOnly = TRUE)[2]

# Load the data
data_info <- read.table(file.path(path, paste0("info_", base_name, ".txt")), col.names = c("INFO"))
print("INFO data loaded")

data_gp <- read.table(file.path(path, paste0("maxGP_", base_name, ".txt")), col.names = c("sample", "MaxGenotypeProbability"))
print("MaxGP data loaded")

# Count the number of SNPs with an INFO value higher than 0.9
count <- sum(data_info$INFO > 0.9)
print(paste("Number of SNPs with INFO value > 0.9:", count))

# Count the number of SNPs with an INFO value higher than 0.95
count <- sum(data_info$INFO > 0.95)
print(paste("Number of SNPs with INFO value > 0.95:", count))

# Count the number of SNPs with an INFO value higher than 0.99
count <- sum(data_info$INFO > 0.99)
print(paste("Number of SNPs with INFO value > 0.99:", count))
### GLOBAL PLOTS ###

# Create the global info and MaxGP histograms
p_info <- ggplot(data_info, aes(x = INFO)) +
  geom_histogram(binwidth = 0.01, color = "black") +
  labs(title = paste("Histogram of INFO (imputation QC) for", base_name), x = "INFO score", y = "Number of SNPs")

p_gp <- ggplot(data_gp, aes(x = MaxGenotypeProbability)) +
  geom_histogram(binwidth = 0.01, color = "black") +
  labs(title = paste("Histogram of maxGP for", base_name), x = "Maximum genotype probability", y = "Number of genotypes")

# Save the plot to a PDF file
ggsave(file.path(path, paste0("info_", base_name, ".pdf")), p_info)
ggsave(file.path(path, paste0("maxGP_", base_name, ".pdf")), p_gp)
