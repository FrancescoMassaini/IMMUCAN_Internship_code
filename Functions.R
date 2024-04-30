filter_common_patients <- function(df1, df2) {
  # Find the common row names (patients) between the two data frames
  common_patients <- intersect(rownames(df1), rownames(df2))
  
  # Filter both data frames to keep only rows corresponding to common patients
  df1_filtered <- df1[common_patients, , drop = FALSE]
  df2_filtered <- df2[common_patients, , drop = FALSE]
  
  # Return a list containing the filtered data frames
  return(list(df1 = df1_filtered, df2 = df2_filtered, common_patients = common_patients))
}

# Remove version number from ENSEMBL ID
remove_gene_version <- function(df_with_gene_version) {
  rownames(df_with_gene_version) <- rownames(df_with_gene_version) %>% 
    str_remove("\\..*$") 
  return(as.data.frame(df_with_gene_version))
}

# Converting from ENSEMBL ID to GENE SYMBOL. Use org.Hs.eg.db db 
EnsemblID_to_GeneSymbol <- function(raw_counts_file){
  entrz <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(raw_counts_file), columns = "SYMBOL", keytype = "ENSEMBL") # keys are the data to overlap. columns is the column to replace and keytype is the column where to find corrispondences   # entrz is a df storing ENSEMBL ID and GENE SYMBOLS. ENSEMBL ID are selected by the ROW NAMES of your raw_counts_file 
  raw_counts_file <- raw_counts_file %>%
    mutate(ENSEMBL = rownames(raw_counts_file)) %>% # add a column
    inner_join(., entrz, by="ENSEMBL") %>% # join, merge df with a common column (ENSEMBL)
    dplyr::filter(!is.na(SYMBOL)) %>%
    distinct(SYMBOL, .keep_all = T) %>% # keep only one variable of the duplicated symbols
    column_to_rownames("SYMBOL") %>% 
    dplyr::select(., -c("ENSEMBL"))
  return(raw_counts_file)  
}

# Raw counts to TPM. Need GENE SYMBOL!
counts_to_TPM <- function(raw_counts_file){
  # log2(TPM+1) 
  TPM <- ADImpute::NormalizeTPM(raw_counts_file, log = T)
  return(as.data.frame(TPM))  
}


# Plots - Heatmap
Compute_Samples_Heatmap <- function(df){
  df <- as.matrix(df)
  
  # Utilizes the NormalizeTPM function from the ADImpute package.
  # Performs normalization on the 'TPM' data frame, considering TPM (Transcripts Per Million) values.
  # Setting 'log = T' indicates that the normalization will be done using the logarithm of TPM values.
  # This is often employed to reduce variance and approximate a normal distribution.
  # Add 1 at the log log(TPM+1) so that 0 values are propertly calculated
  
  sampleDists <- dist(t(df), method = 'euclidean') #Compute distance of the matrix. DO NOT DO IT FOR GENES. dist take into account rows! So to look at the patients you need to transpose. Distance is computed with eucledian metric. Remember you need to transpose (t) the matrix because dist takes rows (and samples are into columns) 
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,    main = "Samples heatmap", sub="", xlab="",
           cex.lab = 1, cex.axis = 1, cex.main = 2, clustering_method = 'complete')
}

# Plots - Dendrogram
Compute_Samples_dendrogram <- function(df){
  sampleTree = hclust(dist(t(df), method = 'euclidean'), method = "complete");
  plot(sampleTree, main = "Samples dendrogram", sub="", xlab="", 
       cex.lab = 1, cex.axis = 1, cex.main = 2)
}


# Correlation
cors <- function(df, cor_type) {
  M <- Hmisc::rcorr(as.matrix(df), type = cor_type)
  Mdf <- map(M, ~data.frame(.x))
  colnames(Mdf[[1]]) <- colnames(M[[1]])
  colnames(Mdf[[2]]) <- colnames(M[[1]])
  colnames(Mdf[[3]]) <- colnames(M[[1]])
  return(Mdf)
}

compute_correlation <- function(deconv_path, deconv_low_value_filter = F, imaging_path, imaging_pattern_to_remove, deconv_pattern_to_remove, significance_level = 0.05, low_corr_threshold = -0.5, high_corr_threshold = 0.5, cor_type = "spearman", imaging_type, deconv_type) {
  
  # Load the data and clean the row names
  imaging <- import(imaging_path)
  rownames(imaging) <- imaging[,1] 
  imaging <- imaging[,-1]
  
  if(imaging_pattern_to_remove != "") rownames(imaging) <- str_remove_all(string = rownames(imaging), pattern = imaging_pattern_to_remove)
  
  deconv <- import(deconv_path)
  rownames(deconv) <- deconv[,1]
  deconv <- deconv[,-1]
  if(deconv_low_value_filter == T) deconv[deconv <= 0.001] <- NA
  
  if(deconv_pattern_to_remove != "") rownames(deconv) <- str_remove_all(string = rownames(deconv), pattern = deconv_pattern_to_remove)
  
  # Find common patients between deconvolution data and the imaging dataset
  common_patients <- intersect(rownames(deconv), rownames(imaging))
  
  # Subset the data for common patients
  imaging_common <- imaging[common_patients,]
  col_clust_imaging <- hclust(dist(na.omit(t(imaging[,-1])), method = "euclidean"))
  deconv_common <- deconv[common_patients, ]
  col_clust_deconv <- hclust(dist(na.omit(t(deconv)), method = "euclidean"))
  
  imaging_common <- imaging_common[, col_clust_imaging$order]
  deconv_common <- deconv_common[, col_clust_deconv$order]
  
  # Merge the datasets
  merged_df <- cbind(imaging_common, deconv_common)
  colnames(merged_df) <- make.unique(colnames(merged_df))
  
  # Compute the correlation
  correlation_matrix <- cors(merged_df, cor_type)
  
  # Clean the correlation matrix based on the significance level
  corr_mat_adjusted <- correlation_matrix
  corr_mat_adjusted[[1]][corr_mat_adjusted[[3]] >= significance_level] <- NA # Put NA to all the r values when p values (which are stored in the 3th position of corr_mat_adjusted) are >= significance_level
  
  # Convert to matrix and apply thresholds for low and high correlations
  corr_mat_adjusted <- as.matrix(corr_mat_adjusted[[1]])
  corr_mat_adjusted[(corr_mat_adjusted <= high_corr_threshold & corr_mat_adjusted >= low_corr_threshold)] <- NA
  
  # Find indices of no NA values in corr_mat_adjusted
  feature1_index <- as.integer(which(!is.na(corr_mat_adjusted), arr.ind = T)[,1])
  feature2_index <- as.integer(which(!is.na(corr_mat_adjusted), arr.ind = T)[,2])
  
  # Take only indices which are not symmetric 
  valid_indices <- feature1_index < feature2_index  
  
  # Map a function to store each correlation value 
  corr_values <- mapply(function(r, c) corr_mat_adjusted[r, c], feature1_index, feature2_index) 
  
  # Create a dataframe to store high correlation pairs
  high_corr_pairs <- data.frame("ID1" = ifelse(rownames(corr_mat_adjusted)[feature1_index] %in% colnames(imaging_common), imaging_type, deconv_type), 
                                "Feature1" = rownames(corr_mat_adjusted)[feature1_index], 
                                "ID2" = ifelse(colnames(corr_mat_adjusted)[feature2_index] %in% colnames(imaging_common), imaging_type, deconv_type), 
                                "Feature2" = colnames(corr_mat_adjusted)[feature2_index],
                                "Corr" = corr_values)
  
  # Filter out high correlation pairs based on valid indices
  high_corr_pairs_filtered <- high_corr_pairs[valid_indices,]
  
  # Splitting correlation based on where they come from
  high_corr_pairs_imaging <- filter(high_corr_pairs_filtered, ID1 == imaging_type, ID2 == imaging_type)
  high_corr_pairs_deconv <- filter(high_corr_pairs_filtered, ID1 == deconv_type, ID2 == deconv_type)
  high_corr_pairs_Mixed <- filter(high_corr_pairs_filtered, (ID1 == imaging_type & ID2 == deconv_type)|(ID1 == deconv_type & ID2 == imaging_type)) %>%
    arrange(desc(Corr))
  
  final_list <- list(high_corr_pairs_Mixed, high_corr_pairs_deconv, high_corr_pairs_imaging, high_corr_pairs_filtered, corr_mat_adjusted, merged_df, correlation_matrix)
  
  return(final_list)
}

#compute_correlation_by_cell_types <- 
#}

plot_correlation_matrix <- function(high_corr_pairs, cor_type, file_name, output_folder, save_png = T, x_label, y_label, png_width = 10, png_height = 10, significance_level = 0.05, low_corr_threshold = -0.5, high_corr_threshold = 0.5) {
  plot <- ggplot(data = high_corr_pairs, aes(x = Feature1, y = Feature2, fill = Corr, label = round(Corr, 2))) +
    geom_tile() +
    labs(fill = paste0(cor_type, "'s\nCorrelation"), title = paste("Correlations", y_label ,"Features vs", x_label),
         subtitle = paste0("Only significant ", cor_type, "'s correlation coefficients shown.\nSignificance level: ", significance_level, ". Low corr threshold: ", low_corr_threshold, ". High corr threshold: ", high_corr_threshold)) +
    scale_fill_gradient2(mid = "#FBFEF9", low = "#0C6291", high = "#A63446", limits = c(-1, 1)) +
    geom_text() +
    theme_classic() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(text = element_text(family = "Roboto"), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab(x_label) + 
    ylab(y_label)
  
  if (save_png == T){
    # Save the plot
    ggsave(filename = paste0(output_folder, file_name, "_" , cor_type, ".png"), plot = last_plot(), width = png_width, height = png_height)
    message(paste("File saved in", file_name))  
  }
  return(plot)
}