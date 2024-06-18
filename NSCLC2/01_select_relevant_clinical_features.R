library(dplyr)
library(stringr)

source("../Functions.R")

clinical = read.csv("../../IMMUCAN_data/NSCLC2/01_Clinical_Data/Daniel/NSCLC2_for_R_update_2024.csv") %>%
  dplyr::select(immucan_id, stage, death, neo_adjuvant_treatment, adjuvant_treatment, simple_histology, long_survivors, TIL_score, side_localization, gender) %>% # select relevat columns for further analysis. 
  # note: time_to_event column is not considered yet here
  column_to_rownames(var = "immucan_id") %>%
  write.csv("../../IMMUCAN_data/NSCLC2/01_Clinical_Data/Daniel/NSCLC2_Daniel_relevant_columns.csv")
