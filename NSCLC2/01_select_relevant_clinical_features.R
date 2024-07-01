library(dplyr)
library(tibble)
library(stringr)

source("../Functions.R")

clinical = read.csv("../../IMMUCAN_data/NSCLC2/01_Clinical_Data/Daniel/NSCLC2_for_R_update_2024.csv") %>%
  dplyr::select(immucan_id, 
                stage, 
                death, 
                # neo_adjuvant_treatment, # interesting but not balanced values
                adjuvant_treatment, 
                simple_histology, # Interesting but not so much balanced and could even be confusing having such different histology types (like adenocarcinoma and squamous cell)
                # long_survivors, # interesting but too many NAs 
                TIL_score, 
                side_localization, 
                gender) %>% # select relevat columns for further analysis.
  dplyr::filter(simple_histology == "Adenocarcinoma") %>% # Take into account only Adenocarcinoma histology 
  dplyr::select(-simple_histology) # delete column after selection

# note: time_to_event column is not considered yet here
clinical = column_to_rownames(clinical, var = "immucan_id")
write.csv(clinical, "../../IMMUCAN_data/NSCLC2/01_Clinical_Data/Daniel/NSCLC2_Daniel_relevant_columns.csv")

