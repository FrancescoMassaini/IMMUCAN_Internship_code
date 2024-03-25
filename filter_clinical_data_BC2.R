library(tidyverse)
library(dplyr)
library(BiocGenerics)
library(stringr)


BC2_clinical <- read.csv("/home/francesco.massaini/Desktop/IMMUCAN_data/BC2/01_Clinical_Data/BC2_colDATA.csv", row.names = 1)

custom_BC2_samples <- read.table("/home/francesco.massaini/Desktop/IMMUCAN_data/BC2/Patient_to_include_from_rna_samples_bc2_to_create_custom_cohort.tsv", header = TRUE,  sep = '	', row.names = NULL)

custom_BC2_samples_unique <- unique(custom_BC2_samples$sample) %>% # Remove duplicates 
  str_replace_all("-",".") # Replace and use the same samples name pattern of all

BC2_clinical <- BC2_clinical[intersect(custom_BC2_samples_unique, rownames(BC2_clinical)),] %>%
  rownames_to_column("Samples")
write.table(BC2_clinical, file = "/home/francesco.massaini/Desktop/IMMUCAN_data/BC2/01_Clinical_Data/BC2_Clinical_data_filtered.txt", quote = F, sep = "\t", row.names = F)

