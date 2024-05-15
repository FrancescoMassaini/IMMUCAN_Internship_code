library(dplyr)
library(tibble)

Daniel_cell_densities = read.csv("/home/francesco.massaini/Desktop/IMMUCAN_data/NSCLC2/01_Imaging_Daniel/IMC_cell_type_densities.csv", header = TRUE)

random_sample = Daniel_cell_densities[sample(nrow(Daniel_cell_densities), 1),"immucan_id"] 
filenames <- list.files("/home/francesco.massaini/Desktop/IMMUCAN_data/NSCLC2/01_Imaging/IMC/single_file_for_patient", pattern=random_sample, full.names=TRUE)

## option 1: sum, every patch
IMC = read.csv(filenames, header = TRUE) 

Daniel_sample = filter(Daniel_cell_densities, immucan_id == random_sample)
Daniel_sample = Daniel_sample[,-2] %>%
  column_to_rownames("immucan_id")
rownames(Daniel_sample) = paste0(rownames(Daniel_sample), "_Daniel")

comparison <- data.frame(as.list(table(IMC$celltypes)))
rownames(comparison) = paste0(random_sample, "_all_patches")

comparison1 = rbind(Daniel_sample, comparison)

## option 2: sum, only tumor
IMC = read.csv(filenames, header = TRUE) %>%
  filter(tumor_patches == TRUE)

Daniel_sample = filter(Daniel_cell_densities, immucan_id == random_sample)
Daniel_sample = Daniel_sample[,-2] %>%
  column_to_rownames("immucan_id")
rownames(Daniel_sample) = paste0(rownames(Daniel_sample), "_Daniel")

comparison <- data.frame(as.list(table(IMC$celltypes)))
rownames(comparison) = paste0(random_sample, "_only_tumor")

comparison2 = rbind(Daniel_sample, comparison)


## option 3: sum, only stroma
IMC = read.csv(filenames, header = TRUE) %>%
  filter(tumor_patches == FALSE) 

Daniel_sample = filter(Daniel_cell_densities, immucan_id == random_sample)
Daniel_sample = Daniel_sample[,-2] %>%
  column_to_rownames("immucan_id")
rownames(Daniel_sample) = paste0(rownames(Daniel_sample), "_Daniel")

comparison <- data.frame(as.list(table(IMC$celltypes)))
rownames(comparison) = paste0(random_sample, "_only_stroma")

comparison3 = rbind(Daniel_sample, comparison)


## option 4: all ROIs average
IMC = read.csv(filenames, header = TRUE) %>%
  group_by(ROI, celltypes) %>%
  summarise(count = n(), .groups = 'drop')

ROI_mean = IMC %>%
  group_by(celltypes) %>%
  summarise(ROI_average = round(mean(count), digits = 0), .groups = "drop")
ROI_mean = as.data.frame(t(ROI_mean))
colnames(ROI_mean) = ROI_mean[1,]
ROI_mean = ROI_mean[-1, , drop=FALSE]  
rownames(ROI_mean) = paste0(random_sample, "_", rownames(ROI_mean)) 
  
Daniel_sample = filter(Daniel_cell_densities, immucan_id == random_sample)
Daniel_sample = Daniel_sample[,-2] %>%
  column_to_rownames("immucan_id")
rownames(Daniel_sample) = paste0(rownames(Daniel_sample), "_Daniel")

comparison4 = rbind(Daniel_sample, ROI_mean)


## option 5: ROI 4 excluded
IMC = read.csv(filenames, header = TRUE) %>%
  filter(ROI != 4) %>%
  group_by(ROI, celltypes) %>%
  summarise(count = n(), .groups = 'drop')

ROI_mean = IMC %>%
  group_by(celltypes) %>%
  summarise(ROI_average_NO_ROI4 = round(mean(count), digits = 0), .groups = "drop")
ROI_mean = as.data.frame(t(ROI_mean))
colnames(ROI_mean) = ROI_mean[1,]
ROI_mean = ROI_mean[-1, , drop=FALSE]  
rownames(ROI_mean) = paste0(random_sample, "_", rownames(ROI_mean)) 

Daniel_sample = filter(Daniel_cell_densities, immucan_id == random_sample)
Daniel_sample = Daniel_sample[,-2] %>%
  column_to_rownames("immucan_id")
rownames(Daniel_sample) = paste0(rownames(Daniel_sample), "_Daniel")

comparison5 = rbind(Daniel_sample, ROI_mean)


df = data.frame("Daniel" = unlist(comparison5[1,]), "Our" = unlist(comparison5[2,]))
df$Our <- as.numeric(df$Our)
df$Daniel <- as.numeric(df$Daniel)
cor(df$Our, df$Daniel)

