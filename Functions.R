filter_common_patients <- function(df1, df2) {
  # Find the common row names (patients) between the two data frames
  common_patients <- intersect(rownames(df1), rownames(df2))
  
  # Filter both data frames to keep only rows corresponding to common patients
  df1_filtered <- df1[common_patients, , drop = FALSE]
  df2_filtered <- df2[common_patients, , drop = FALSE]
  
  # Return a list containing the filtered data frames
  return(list(df1 = df1_filtered, df2 = df2_filtered, common_patients = common_patients))
}