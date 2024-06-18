library(dplyr)
library(stringr)

source("../Functions.R")

input_folder <- "../../IMMUCAN_data/NSCLC2/01_Raw_Counts/star_counts"
txt_input_files <- list.files(path = input_folder, pattern = "\\.txt$", full.names = TRUE)

# Storing sample names from file names 
files_names <- gsub(".txt$", "", basename(txt_input_files))
samples_names = str_split_i(files_names, "_", 1)

# Merging all the patients files and removing gene version
merged_counts <- read.delim(txt_input_files[1], header = 0, row.names = 1)
colnames(merged_counts) <- samples_names[1]

for (i in 2:length(txt_input_files)) {
  counts <- read.delim(txt_input_files[i], header = 0, row.names = 1)
  colnames(counts) <- samples_names[i]
  # Unisci le colonne del nuovo data frame a merged_counts
  merged_counts <- cbind(merged_counts, counts)
}

cat("Sum of counts across samples: \n")
colSums(merged_counts) # Not yet normalized to TPM

# Deconvolution methods already take into account low expressed and low variance genes so we can directly compute TPM
TPM_forDeconv = merged_counts %>%  
  filter(!rownames(.) %in% c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")) %>%
  remove_gene_version() %>%
  EnsemblID_to_GeneSymbol() %>% # Deconv requires gene symbols
  counts_to_TPM(., log=T)

write.table(TPM_forDeconv, "../../IMMUCAN_data/NSCLC2/02_TPM/NSCLC2_TPM_forDecov.txt", sep = "\t")

# For Differential Expression analysis let's keep ENSEMBL ID and counts
merged_counts_forDE = merged_counts %>%  
  filter(!rownames(.) %in% c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")) %>%
  remove_gene_version() 

write.table(merged_counts_forDE, "../../IMMUCAN_data/NSCLC2/01_Raw_Counts/NSCLC2_all_patients_counts_forDE.txt", sep = "\t")
