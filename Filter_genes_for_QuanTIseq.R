# Deleting (not)duplicated file to avoid problems with QuanTIseq
library(tidyverse)
library(dplyr)

# path <- c("/home/francesco.massaini/Desktop/IMMUCAN_data/BC1/02_TPM/BC1_TPM_Mounim.txt",
#                "/home/francesco.massaini/Desktop/IMMUCAN_data/BC2/02_TPM/BC2_TPM_Mounim.txt", 
#                "/home/francesco.massaini/Desktop/IMMUCAN_data/BC3/02_TPM/BC3_TPM_Mounim.txt",
#                "/home/francesco.massaini/Desktop/IMMUCAN_data/NSCLC/02_TPM/NSCLC_TPM_Mounim.txt")
path <- ("/home/francesco.massaini/Desktop/IMMUCAN_data/BC2/02_TPM/BC2_TPM_fromRdata.txt")
          
files_names <- gsub(".txt$", "", basename(path))
folders_names <- strsplit(files_names, "_")

for (i in 1:length(path)) {
  TPM <- read.table(path[i], header = TRUE,  sep = '	', row.names = 1)
  TPM <- TPM[rownames(TPM)%in%rownames(immunedeconv::dataset_racle$expr_mat),] %>%
    rownames_to_column("Genes") 
  write.table(TPM, file = paste0("/home/francesco.massaini/Desktop/IMMUCAN_data/", folders_names[[i]][1], "/02_TPM/", files_names[i],"_genes_filtered_for_Quantiseq.txt"), quote = F, sep = "\t", row.names = F)
}