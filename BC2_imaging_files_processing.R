library(stringr)
library(dplyr)

IF_path = "../IMMUCAN_data/BC2/01_Imaging/IF/"
IF_dirs = c(paste0(IF_path, "IF3/")) ####
# IF_dirs = c(paste0(IF_path, "IF1/"), paste0(IF_path, "IF2/"), paste0(IF_path, "IF3/"))) IT CRASHES DUE TO A LACK OF RAM

all = data.frame()
for (IF_dir in IF_dirs) {
  print(paste("Processing directory: ", IF_dir))
  
  all_files <- list.files(IF_dir)
  file_info <- file.info(file.path(IF_dir, all_files))
  
  # Filter only files (not directories)
  only_files <- all_files[!file_info$isdir]
  
  # Right phenotype key file based on the current IF panel
  phenotype_keys = read.csv(paste0("../IMMUCAN_data/Phenotype_keys_from_GitHub/phenotype_key_", str_extract(string = IF_dir, pattern = "IF\\d"),".csv"))
  print(paste("Phenotype keys: ", paste0("../IMMUCAN_data/Phenotype_keys_from_GitHub/phenotype_key_", str_extract(string = IF_dir, pattern = "IF\\d"),".csv")))
  
  for (file in only_files) {
    print(paste("Processing: ", file))
    if (file.info(paste0(IF_dir, file))$size > 0) {
      df <- read.delim(paste0(IF_dir, file), sep = "\t") %>%
        select("cell.ID", "nucleus.x", "nucleus.y", "tissue.type", "phenotype") %>%
        left_join(phenotype_keys, by = "phenotype") %>%
        mutate(celltype = coalesce(celltype, phenotype))
      colnames(df) = str_replace_all(string = colnames(df), pattern = "\\.", replacement = "_")
      df["file_name"] <- file
      df["sample_panel_ID"] <- str_split_i(string = file, pattern = "_#", i = 1) #### DA CORREGGERE SOTTO
      df["sample_ID"] <- str_split_i(string = file, pattern = "-IF\\d", i = 1)
      df["patient_ID"] <- str_split_i(string = file, pattern = "-FIXT", i = 1)
      df["panel_ID"] <- str_extract(string = file, pattern = "IF\\d-0\\d") ##### DA CORREGGERE SOTTO
      df["panel"] <- str_extract(string = file, pattern = "IF\\d")
      all <- rbind(all, df)
    } else {
      print(paste("File is empty, skipping:", file))
    }
  }
  write.csv(all, paste0(IF_dir, str_extract(string = IF_dir, pattern = "IF\\d"), "_all_patients_cell_coord.csv"), row.names = FALSE)
}

# Optional: Filter patients from a list
## BC2 - From Andrea BC drive
all = read.csv("../IMMUCAN_data/BC2/01_Imaging/IF/IF2/IF2_all_patients_cell_coord.csv")
IF_dir = "../IMMUCAN_data/BC2/01_Imaging/IF/IF2/"
current_IF = "IF2"

tumor_IF_all = all %>%
  filter(tissue_type == "tumor")
tumor_celltype_freq = setNames(as.data.frame(table(tumor_IF_all[,"sample_panel_ID"], tumor_IF_all[,"celltype"])), c("sample_panel_ID", "celltypes", "Freq"))  %>%
  mutate(tissue_type = "tumor") ### DA CORREGGERE SOTTO

write.csv(tumor_celltype_freq, paste0(IF_dir, current_IF,"_all_patients_tumor_cell_freq.csv"))

stroma_IF_all = all %>%
  filter(tissue_type == "stroma")
stroma_celltype_freq = setNames(as.data.frame(table(stroma_IF_all[,"sample_panel_ID"], stroma_IF_all[,"celltype"])), c("sample_panel_ID", "celltypes", "Freq")) %>%
  mutate(tissue_type = "stroma")
write.csv(stroma_celltype_freq, paste0(IF_dir, current_IF, "_all_patients_stroma_cell_freq.csv"))



cat("Total rows before filtering ", nrow(all), "\n", "Total samples", length(unique(all$sample_ID)))
IF1_patients_list = readxl::read_excel("../IMMUCAN_data/BC2/01_Reference_files_from_Andrea_BC2_TNBC/List_samples_data_type/mif1_samples_bc2_custom_cohort.xlsx") %>%
  mutate(panel = "IF1") %>%
  rename(tissue_type_Andrea = tissue_type)
cat("IF1 numb of patients ", nrow(IF1_patients_list))
IF2_patients_list = readxl::read_excel("../IMMUCAN_data/BC2/01_Reference_files_from_Andrea_BC2_TNBC/List_samples_data_type/mif2_samples_bc2_custom_cohort.xlsx") %>%
  mutate(panel = "IF2") %>%
  rename(tissue_type_Andrea = tissue_type)
cat("IF2 numb of patients ", nrow(IF2_patients_list))
IF3_patients_list = readxl::read_excel("../IMMUCAN_data/BC2/01_Reference_files_from_Andrea_BC2_TNBC/List_samples_data_type/mif3_samples_bc2_custom_cohort.xlsx") %>%
  mutate(panel = "IF3") %>%
  rename(tissue_type_Andrea = tissue_type)
cat("IF3 numb of patients ", nrow(IF3_patients_list))
IFs = rbind(IF1_patients_list, IF2_patients_list, IF3_patients_list)
cat("all IFs total rows ", nrow(IF1_patients_list))
# write.csv(IFs, "../IMMUCAN_data/BC2/01_Reference_files_from_Andrea_BC2_TNBC/List_samples_data_type/ALL_mIFs_samples_bc2_custom_cohort.xlsx") ####
####

IF_patients_filtered = all %>%
  inner_join(IFs, by = c("sample_ID" = "sample", "panel" = "panel"))
cat("Total rows after filtering ", nrow(IF_patients_filtered), "\n", "Total samples", length(unique(IF_patients_filtered$sample_ID)))
cat("Diff:", setdiff(IF3_patients_list$sample, unique(IF_patients_filtered$sample_ID)))
#write.csv(IF_patients_filtered, "../IMMUCAN_data/BC2/01_Imaging/IF/IF1/IF1_filtered_patients_Andrea.csv") ####
write.csv(IF_patients_filtered, "../IMMUCAN_data/BC2/01_Imaging/IF/IF2/IF2_filtered_patients_Andrea.csv") ####
# write.csv(IF_patients_filtered, "../IMMUCAN_data/BC2/01_Imaging/IF/IF3/IF3_filtered_patients_Andrea.csv") ####


tumor_IF_filtered = IF_patients_filtered %>%
  filter(tissue_type == "tumor")
tumor_celltype_freq = setNames(as.data.frame(table(tumor_IF_filtered[,"sample_panel_ID"], tumor_IF_filtered[,"celltype"])), c("sample_panel_ID", "celltypes", "Freq"))
write.csv(tumor_celltype_freq, paste0(IF_dir, current_IF, "_filtered_patients_Andrea_tumor_cell_freq.csv"))

stroma_IF_filtered = IF_patients_filtered %>%
  filter(tissue_type == "stroma")
stroma_celltype_freq = setNames(as.data.frame(table(stroma_IF_filtered[,"sample_panel_ID"], stroma_IF_filtered[,"celltype"])), c("sample_panel_ID", "celltypes", "Freq"))
write.csv(stroma_celltype_freq, paste0(IF_dir, current_IF, "_filtered_patients_Andrea_stroma_cell_freq.csv"))













#### Parte nuova
library(stringr)
library(dplyr)
IF_path = "../IMMUCAN_data/BC2/01_Imaging/IF/"
IF_dirs = c(paste0(IF_path, "IF1/"), paste0(IF_path, "IF2/"), paste0(IF_path, "IF3/"))
IF1_patients_list = readxl::read_excel("../IMMUCAN_data/BC2/01_Reference_files_from_Andrea_BC2_TNBC/List_samples_data_type/mif1_samples_bc2_custom_cohort.xlsx") %>%
  mutate(panel = "IF1")
cat("IF1 numb of patients ", nrow(IF1_patients_list))
IF2_patients_list = readxl::read_excel("../IMMUCAN_data/BC2/01_Reference_files_from_Andrea_BC2_TNBC/List_samples_data_type/mif2_samples_bc2_custom_cohort.xlsx") %>%
  mutate(panel = "IF2")
cat("IF2 numb of patients ", nrow(IF2_patients_list))
IF3_patients_list = readxl::read_excel("../IMMUCAN_data/BC2/01_Reference_files_from_Andrea_BC2_TNBC/List_samples_data_type/mif3_samples_bc2_custom_cohort.xlsx") %>%
  mutate(panel = "IF3")
cat("IF3 numb of patients ", nrow(IF3_patients_list))
IFs = rbind(IF1_patients_list, IF2_patients_list, IF3_patients_list)
cat("all IFs total rows ", nrow(IF1_patients_list))
write.csv(IFs, "../IMMUCAN_data/BC2/01_Reference_files_from_Andrea_BC2_TNBC/List_samples_data_type/ALL_mIFs_samples_bc2_custom_cohort.xlsx") ####
all = data.frame()
for (IF_dir in IF_dirs) {
  print(paste("Processing directory: ", IF_dir))
  current_IF = str_extract(string = IF_dir, pattern = "IF\\d")
  all_files <- list.files(IF_dir)
  file_info <- file.info(file.path(IF_dir, all_files))
  # Filter only files (not directories)
  only_files <- all_files[!file_info$isdir]
  # Right phenotype key file based on the current IF panel
  phenotype_keys = read.csv(paste0("../IMMUCAN_data/Phenotype_keys_from_GitHub/phenotype_key_", current_IF,".csv"))
  print(paste("Phenotype keys: ", paste0("../IMMUCAN_data/Phenotype_keys_from_GitHub/phenotype_key_", current_IF,".csv")))
  for (file in only_files) {
    print(paste("Processing: ", file))
    if (file.info(paste0(IF_dir, file))$size > 0) {
      df <- read.delim(paste0(IF_dir, file), sep = "\t") %>%
        select("cell.ID", "nucleus.x", "nucleus.y", "tissue.type", "phenotype") %>%
        left_join(phenotype_keys, by = "phenotype")
      colnames(df) = str_replace_all(string = colnames(df), pattern = "\\.", replacement = "_")
      df["file_name"] <- file
      df["sample_panel_ID"] <- str_split_i(string = file, pattern = "#", i = 1)
      df["sample_ID"] <- str_split_i(string = file, pattern = "-IF\\d", i = 1)
      df["patient_ID"] <- str_split_i(string = file, pattern = "-FIXT", i = 1)
      df["panel_ID"] <- str_extract(string = file, pattern = "IF\\d")
      df["panel"] <- str_extract(string = file, pattern = "IF\\d")
      all <- rbind(all, df)
    } else {
      print(paste("File is empty, skipping:", file))
    }
  }
  write.csv(all, paste0(IF_dir, "_all_patients_cell_coord.csv"), row.names = FALSE)
  tumor_IF_all = all %>%
    filter(tissue_type == "tumor")
  tumor_celltype_freq = as.data.frame(table(tumor_IF_all$sample_panel_ID, tumor_IF_all$celltypes))
  
  write.csv(tumor_celltype_freq, paste0(IF_dir, "_all_patients_tumor_cell_freq.csv"))
  
  stroma_IF_all = all %>%
    filter(tissue_type == "stroma")
  stroma_celltype_freq = as.data.frame(table(stroma_IF_all$sample_panel_ID, stroma_IF_all$celltypes))
  write.csv(stroma_celltype_freq, paste0(IF_dir, "_all_patients_stroma_cell_freq.csv"))
  
  # Optional: Filter patients from a list
  ## BC2 - From Andrea BC drive
  cat("Total rows before filtering ", nrow(all), "\n", "Total samples", length(unique(all$sample_ID)))
  ####
  IF_patients_filtered = all %>%
    inner_join(IFs, by = c("sample_ID" = "sample", "panel" = "panel"))
  cat("Total rows after filtering ", nrow(IF_patients_filtered), "\n", "Total samples", length(unique(IF_patients_filtered$sample_ID)))
  ifelse(test = current_IF == "IF1",
         cat("Diff:", setdiff(IF1_patients_list$sample, unique(IF_patients_filtered$sample_ID))),
         ifelse(test = current_IF == "IF2",
                cat("Diff:", setdiff(IF2_patients_list$sample, unique(IF_patients_filtered$sample_ID))),
                cat("Diff:", setdiff(IF3_patients_list$sample, unique(IF_patients_filtered$sample_ID))),
         ))
  write.csv(IF_patients_filtered, paste0(IF_dir, "_filtered_patients_Andrea_cell_coord.csv")) ####
  tumor_IF_filtered = IF_patients_filtered %>%
    filter(tissue_type == "tumor")
  tumor_celltype_freq = as.data.frame(table(tumor_IF_filtered$sample_panel_ID, tumor_IF_filtered$celltypes))
  write.csv(IF_patients_filtered, paste0(IF_dir, "_filtered_patients_Andrea_tumor_cell_freq.csv"))
  stroma_IF_filtered = IF_patients_filtered %>%
    filter(tissue_type == "stroma")
  stroma_celltype_freq = as.data.frame(table(stroma_IF_filtered$sample_panel_ID, stroma_IF_filtered$celltypes))
  write.csv(IF_patients_filtered, paste0(IF_dir, "_filtered_patients_Andrea_stroma_cell_freq.csv"))
}