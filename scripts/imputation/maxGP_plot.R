# This R script generates a pdf plot that representes the maximum GP distribution of a TXT file generated from a VCF with the command: 
# bcftools query -f '[%GP\n]' medcov_autosomes.vcf | awk -F',' '{print ($1>$2)?($1>$3?$1:$3):($2>$3?$2:$3)}'

# The input file is the TXT file table with the number of SNPs with each frequency of missing data generated with bcftools query. 

# Usage: Rscript maxGP_plot.R <input_path> <input_basename>
library(ggplot2)

# Get the path and base name from the command-line arguments
path <- commandArgs(trailingOnly = TRUE)[1]
base_name <- commandArgs(trailingOnly = TRUE)[2]

# Load the data
data_gp <- read.table(file.path(path, paste0(base_name, ".txt")), col.names = c("Sample", "MaxGenotypeProb"))
print("MaxGP data loaded")

# Create the global maxGP histogram
p_gp <- ggplot(data_gp, aes(x = MaxGenotypeProb)) +
  geom_histogram(binwidth = 0.01, color = "black") +
  labs(title = paste("Histogram of maximum genotype probabilities for", base_name), x = "Maximum genotype probability", y = "Number of SNPs")

# Save the plot to a PDF file
ggsave(file.path(path, paste0(base_name, ".pdf")), p_gp)
print("Global maxGP histogram saved. Calculating individual plots...")


# Get the unique samples
samples <- unique(data_gp$Sample)

# Calculate the number of plot groups
num_groups <- ceiling(length(samples) / 4)

# Loop over the plot groups
for (i in seq_len(num_groups)) {
  # Get the samples for the current group
  group_samples <- samples[((i - 1) * 4 + 1):min(i * 4, length(samples))]
  
  # Subset the data for the current group
  data_group <- subset(data_gp, Sample %in% group_samples)
  
  # Create the histogram for the current group
  p_gp_group <- ggplot(data_group, aes(x = MaxGenotypeProb)) +
    geom_histogram(binwidth = 0.01, color = "black") +
    facet_wrap(~ Sample, ncol = 2) +
    labs(title = paste("Histogram of maximum genotype probabilities for", base_name, "group", i), x = "Maximum genotype probability", y = "Number of SNPs")
  
  # Save the plot to a PDF file
  ggsave(file.path(path, paste0(base_name, "_group", i, ".pdf")), p_gp_group)
}
