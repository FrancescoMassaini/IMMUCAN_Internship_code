library(stringr)
library(dplyr)

phenotype_keys = read.csv("../IMMUCAN_data/Phenotype_keys_from_GitHub/phenotype_key_IF2.csv")
files <- list.files("../IMMUCAN_data/BC2/IF2/")
all = data.frame()
for (file in files) {
  df = read.delim(paste0("../IMMUCAN_data/BC2/IF2/", file), sep = "\t") %>%
    select("cell.ID", "nucleus.x", "nucleus.y", "tissue.type", "phenotype") %>%
    left_join(phenotype_keys, by = "phenotype")
  #df["file.name"] = file
  #df["sample.ID"] = str_split_i(string = file, pattern = "_#_", i = 1)
  #df["patient.ID"] = str_split_i(string = file, pattern = "-FIXT", i = 1)
  #df["panel"] = str_extract(string = file, pattern = "IF\\d")
  #all = rbind(all, df)
}
