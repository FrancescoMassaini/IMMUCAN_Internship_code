library(tidyverse)
library(dplyr)
library(BiocGenerics)
library(stringr)

all_BCs <- read.table("/home/francesco.massaini/Desktop/IMMUCAN_data/all_BCs_merged_BC-20032024_tpm.txt", header = TRUE,  sep = '	', row.names = 1)
NSCLC <- read.table("/home/francesco.massaini/Desktop/IMMUCAN_data/NSCLC/NSCLC_tpm_Mounim.txt", header = TRUE,  sep = '	', row.names = 1) %>%
  rownames_to_column("Genes") 

BC1_samples <- grep(pattern = "BC1", colnames(all_BCs), value = T)
BC3_samples <- grep(pattern = "BC3", colnames(all_BCs), value = T)

# BC2 needs to be filtered with custom patients (specified in Patient_to_include_from_rna_samples_bc2_to_create_custom_cohort.tsv) 
custom_BC2_samples <- read.table("/home/francesco.massaini/Desktop/IMMUCAN_data/BC2/Patient_to_include_from_rna_samples_bc2_to_create_custom_cohort.tsv", header = TRUE,  sep = '	', row.names = NULL)
custom_BC2_samples_unique <- unique(custom_BC2_samples$sample) %>% # Remove duplicates 
  str_replace_all("-",".") %>%
  paste0(.,".RNA.01") # Replace and use the same samples name pattern of all
BC2_samples <- grep(pattern = "BC2", colnames(all_BCs), value = T) %>%
  intersect(., custom_BC2_samples_unique)

# Save distinct files (with filtered BC2) 
BC1 <- all_BCs[BC1_samples] %>%
  rownames_to_column("Genes")
write.table(BC1, file = "/home/francesco.massaini/Desktop/IMMUCAN_data/BC1/02_TPM/BC1_TPM_Mounim.txt", quote = F, sep = "\t", row.names = F)
BC2 <- all_BCs[BC2_samples] %>%
  rownames_to_column("Genes") 
write.table(BC2, file = "/home/francesco.massaini/Desktop/IMMUCAN_data/BC2/02_TPM/BC2_TPM_Mounim.txt", quote = F, sep = "\t", row.names = F)
BC3 <- all_BCs[BC3_samples] %>%
  rownames_to_column("Genes") 
write.table(BC3, file = "/home/francesco.massaini/Desktop/IMMUCAN_data/BC3/02_TPM/BC3_TPM_Mounim.txt", quote = F, sep = "\t", row.names = F)

write.table(NSCLC, file = "/home/francesco.massaini/Desktop/IMMUCAN_data/NSCLC/02_TPM/NSCLC_TPM_Mounim.txt", quote = F, sep = "\t", row.names = F)

  
