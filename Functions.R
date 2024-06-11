library(ggplot2)
library(ggdendro)
library(dendextend)

dfs_same_patients_same_order <- function(expr_df, sample_feat_df, expr_patients_on_rows = FALSE, clinical_patients_on_rows = TRUE) {
  # Initial transpose of the expression dataframe
  
  df_transposed <- as.data.frame(t(expr_df))

  # Merging with patient metadata
  df_merged <- merge(df_transposed, sample_feat_df, by = 0) %>% 
    column_to_rownames(var = "Row.names")
  
  # Select columns excluding metadata for the expression dataframe
  df_expr <- dplyr::select(df_merged, -colnames(sample_feat_df))
  
  # Select only metadata columns for the clinical dataframe
  df_clinical <- dplyr::select(df_merged, colnames(sample_feat_df))
  
  # Final transpose of dataframes based on parameters
  if (expr_patients_on_rows == FALSE) {
    df_expr <- as.data.frame(t(df_expr)) 
  }
  
  if (clinical_patients_on_rows == FALSE) {
    df_clinical <- as.data.frame(t(df_clinical))
  }
  
  # Return both dataframes
  list(expression_data = df_expr, clinical_data = df_clinical)
}


# filter_common_patients <- function(df1, df2) {
#   # Find the common row names (patients) between the two data frames
#   common_patients <- intersect(rownames(df1), rownames(df2))
#   
#   # Filter both data frames to keep only rows corresponding to common patients
#   df1_filtered <- df1[common_patients, , drop = FALSE]
#   df2_filtered <- df2[common_patients, , drop = FALSE]
#   
#   # Return a list containing the filtered data frames
#   return(list(df1 = df1_filtered, df2 = df2_filtered, common_patients = common_patients))
# }

 
# # Remove version number from ENSEMBL ID
remove_gene_version <- function(df_with_gene_version) {
  rownames(df_with_gene_version) <- rownames(df_with_gene_version) %>%
    str_remove("\\..*$")
  return(as.data.frame(df_with_gene_version))
}

# Converting from ENSEMBL ID to GENE SYMBOL. Use org.Hs.eg.db db 
EnsemblID_to_GeneSymbol <- function(raw_counts_file, aggregation_method = "mean"){

  entrz <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(raw_counts_file), columns = "SYMBOL", keytype = "ENSEMBL") # keys are the data to overlap. columns is the column to replace and keytype is the column where to find corrispondences   # entrz is a df storing ENSEMBL ID and GENE SYMBOLS. ENSEMBL ID are selected by the ROW NAMES of your raw_counts_file 
  raw_counts_file_gene_symbol <- raw_counts_file %>%
    mutate(ENSEMBL = rownames(raw_counts_file)) %>% # add a column
    inner_join(., entrz, by="ENSEMBL") %>% # join, merge df with a common column (ENSEMBL)
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::select(-ENSEMBL) %>%
    group_by(SYMBOL)
  
  # Duplicated gene symbols aggregation based on the specified method 
  if (aggregation_method == "mean") {
    raw_counts_file_gene_symbol <- raw_counts_file_gene_symbol %>%
      summarise(across(everything(), mean), .groups = 'drop') %>% # take the average of duplicated rows
      column_to_rownames("SYMBOL")
  } else if (aggregation_method == "sum") {
    raw_counts_file_gene_symbol <- raw_counts_file_gene_symbol %>%
      summarise(across(everything(), sum), .groups = 'drop') %>% # take the sum of duplicated rows
      column_to_rownames("SYMBOL")
  }
  
  colnames(raw_counts_file_gene_symbol) = gsub(pattern = "_1$", replacement = "", x = colnames(raw_counts_file_gene_symbol))
  
  cat("Number of genes before gene symbol conversion:", nrow(raw_counts_file), "\n", "Number of genes after conversion:", nrow(raw_counts_file_gene_symbol), "\n", "Difference:", nrow(raw_counts_file) - nrow(raw_counts_file_gene_symbol))
  cat("\n")
  return(raw_counts_file_gene_symbol)  
}

# Raw counts to TPM. Need GENE SYMBOL!
counts_to_TPM <- function(raw_counts_file, log = T){
  # log2(TPM+1) 
  TPM <- ADImpute::NormalizeTPM(raw_counts_file, log = log)
  return(as.data.frame(TPM))  
}

# Define a function to filter genes based on average expression and zero count frequency
filter_low_expr_genes <- function(data, avg_expr_threshold = 10, zero_count_percent_threshold = 0.9) {
  # Calculate the mean expression for each gene
  # Calculate the frequency of zero counts for each gene
  data %>%
    mutate(
      mean_expression = rowMeans(.),  # Add a column for the mean expression of each gene
      zero_count_frequency = rowSums(. == 0) / ncol(.)  # Add a column for the frequency of zero counts
    ) %>%
    filter(
      mean_expression >= avg_expr_threshold,  # Filter genes with mean expression below the threshold
      zero_count_frequency <= zero_count_percent_threshold  # Filter genes with a high frequency of zero counts
    ) %>%
    dplyr::select(-mean_expression, -zero_count_frequency)  # Select only the original count columns for the output
}

# Define a function to filter genes based on variance
filter_low_variance_genes <- function(data, quantile_threshold = 0.05) {
  
  # Calculate the variance for each gene and filter
  df = as.data.frame(rowVars(as.matrix(data)))
  colnames(df) = "variance"
  df = rownames_to_column(df, var = "genes")
  
  # Genes to keep
  keep = df[df$variance > quantile(df$variance, prob = quantile_threshold), "genes"]
    
  data = data[keep, ]
  return(data)
}

# Plots - Distribution
compute_distribution <- function(data, plot_title, xlab, ylab = "Density", use_log = TRUE) {
  # Convert row names to a column
  data_with_genes <- data %>%
    rownames_to_column(var = "Genes")
  
  # Pivot data to long format
  data_long <- pivot_longer(data_with_genes, 
                            cols = -Genes,
                            names_to = "Patient",
                            values_to = "Counts")
  
  # Calculate log2(value+1) if use_log is TRUE, otherwise use the direct value
  if (use_log) {
    # Create density plot
    plot <- ggplot(data_long, aes(x=log2(Counts+1))) +
      geom_density() +
      theme_bw() +
      labs(title = plot_title, 
           x = xlab, 
           y = ylab)
  }
  else {
    plot <- ggplot(data_long, aes(x=Counts)) +
      geom_density() +
      theme_bw() +
      labs(title = plot_title, 
           x = xlab, 
           y = ylab)
  }
  
  return(plot)
}

# Plots - Heatmap
## patients need to start on columns
Compute_Samples_Heatmap <- function(df, sample_feat_df = NULL, main_title = NULL){
  # Ensure same patients order
  tmp = dfs_same_patients_same_order(expr_df = df, sample_feat_df = sample_feat_df, expr_patients_on_rows = TRUE, clinical_patients_on_rows = TRUE)
  df = tmp[["expression_data"]]
  sample_feat_df = tmp[["clinical_data"]]
  
  # Compute Heatmap
  sampleDists <- dist(df, method = 'euclidean') #Compute distance of the matrix. DO NOT DO IT FOR GENES. dist take into account rows! So to look at the patients you need to transpose. Distance is computed with eucledian metric. Remember you need to transpose (t) the matrix because dist takes rows (and samples are into columns) 
  sampleDistMatrix <- as.matrix(sampleDists)
  
  pheatmap(sampleDistMatrix, main = main_title,
           annotation_col = sample_feat_df,
           color = hcl.colors(50, "BluYl"),
           show_colnames = T, show_rownames = T,
           fontsize_col = 5, fontsize_row = 5, fontsize = 10, clustering_method = 'complete',
           legend_breaks = c(min(sampleDistMatrix), max(sampleDistMatrix)),
           legend_labels = c("Similar", "Distant"),
           border = NA)
}


# Plots - Dendrogram
Compute_Samples_dendrogram <- function(df, sample_feat_df, feat_column = NULL, main_title = NULL, palette_colors = "Dark2"){
  # Ensure same patients order
  tmp = dfs_same_patients_same_order(expr_df = df, sample_feat_df = sample_feat_df, expr_patients_on_rows = TRUE, clinical_patients_on_rows = TRUE)
  df = tmp[["expression_data"]]
  sample_feat_df = tmp[["clinical_data"]]
  
  # Compute Dendrogram
  sampleTree = hclust(dist(df, method = 'euclidean'), method = "complete")
  
  dend_expr <- as.dendrogram(sampleTree)
  tree_labels<- dendro_data(dend_expr, type = "rectangle")
  tree_labels$labels <- cbind(tree_labels$labels, feat = as.factor(sample_feat_df[,feat_column]))
  
  ggplot() +
    geom_segment(data = segment(tree_labels), aes(x=x, y=y, xend=xend, yend=yend))+
    geom_segment(data = tree_labels$segments %>%
                 filter(yend == 0) %>%
                 left_join(tree_labels$labels, by = "x"), aes(x=x, y=y.x, xend=xend, yend=yend, color = feat)) +
    # geom_text(data = label(tree_labels), aes(x=x, y=y, label=label, colour = feat, hjust=0), size=1) +
    coord_flip() +
    scale_y_reverse(expand=c(0.2, 0)) +
    scale_colour_brewer(palette = palette_colors, name = feat_column) + 
    theme_dendro()
}

# Plots - Violin plot
Wrapped_violin_plot_by_clinical_feature <- function(df_long, metadata_df, metadata_var, facet_wrap_var){
  metadata_df = rownames_to_column(metadata_df, var = "Sample")
  tmp_clinical = metadata_df[, c("Sample", metadata_var)]
  colnames(tmp_clinical) = c("Sample", "clin")
  merged <- left_join(df_long, tmp_clinical, by = "Sample")
  
  facet_formula <- as.formula(paste("~", facet_wrap_var))
  
  p <- ggplot(merged, aes(x = clin, y = Value, fill = clin)) +
    geom_violin(trim = FALSE) +
    geom_jitter(width = 0.1, alpha = 0.2, size = 1.5, color = "black") +
    labs(title = paste("Violin plot of", clin_var), y = "IC Level", x = clin_var, fill = clin_var) +
    facet_wrap(facet_formula) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
}

# Plots - Only significant violin plot
significant_violin_plot = function(df_long, metadata_df, metadata_var, p_threshold = 0.05, main_title){
  metadata_df = rownames_to_column(metadata_df, var = "Sample")
  tmp_clinical = metadata_df[, c("Sample", metadata_var), drop = FALSE]
  colnames(tmp_clinical) = c("Sample", "clin")
  merged <- left_join(df_long, tmp_clinical, by = "Sample")
  
  test_result <- aov(Value ~ clin, data = merged)
  summary_test_result <- summary(test_result)
  p_value = summary_test_result[[1]]$'Pr(>F)'[1]
  
  if (!is.na(p_value) && p_value < p_threshold) {  # Check p-value
    p <- ggplot(merged, aes(x = clin, y = Value, fill = clin)) +
      geom_violin(trim = FALSE) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, color = "black") +
      labs(title = paste(main_title, "p value: ", round(p_value, digits = 2)), y = "IC Level", x = clin_var, fill = clin_var) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)
  }
  else{
    cat("No significant difference for: ", metadata_var, "p value:", p_value ,"\n")
  }
}
  
# Plots - Correlation plot from two different dfs
## dfs must have the same rownames and features on columns 
create_correlation_plot <- function(data1, data2, correlation_type = "pearson", r_threshold = 0.7, p_threshold = 0.05) {
  # Merge the two datasets
  merged_data = cbind(data1, data2)
  
  # Compute correlation matrix using Hmisc package
  cor_matrix = Hmisc::rcorr(as.matrix(merged_data), type = correlation_type)
  
  # Extract correlation coefficients and p-values
  correlation_data = cor_matrix$r %>%
    as.data.frame() %>%
    rownames_to_column(var = "variable1") %>%
    pivot_longer(cols = -variable1, names_to = "variable2", values_to = "correlation") %>%
    filter(variable1 %in% colnames(data1) & variable2 %in% colnames(data2))
  
  p_value_data = cor_matrix$P %>%
    as.data.frame() %>%
    rownames_to_column(var = "variable1") %>%
    pivot_longer(cols = -variable1, names_to = "variable2", values_to = "p_value") %>%
    filter(variable1 %in% colnames(data1) & variable2 %in% colnames(data2))
  
  # Join the correlation and p-value data
  correlation_analysis = left_join(correlation_data, p_value_data, by = c("variable1", "variable2")) %>%
    filter((correlation < -r_threshold | correlation > r_threshold) & !is.na(p_value) & p_value <= p_threshold)
  
  # Plot the results using ggplot2
  plot = ggplot(data = correlation_analysis, aes(x = variable1, y = variable2, fill = correlation, label = round(correlation, digits = 2))) +
    geom_tile() +
    geom_text(color = "black", size = 3, vjust = 1.5) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(plot)
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