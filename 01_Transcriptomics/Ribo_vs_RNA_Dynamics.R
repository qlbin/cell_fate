# ==============================================================================
# Script: Ribo_vs_RNA_Dynamics.R
# Purpose: Merge Ribo-seq and RNA-seq matrices and calculate Translational Efficiency (TE)
# ==============================================================================
source("00_utils.R")
suppressPackageStartupMessages({
  library(data.table)
  library(MetBrewer)
})

# ------------------------------------------------------------------------------
# 1. Batch Load & Merge Matrix
# ------------------------------------------------------------------------------
merge_feature_files <- function(folder_path, file_pattern = "\\.txt\\.gz$") {
  
  file_list <- list.files(path = folder_path, pattern = file_pattern, full.names = TRUE)
  if(length(file_list) == 0) stop("No files found in the specified directory.")
  
  read_and_format <- function(file) {
    dt <- fread(file, check.names = TRUE)
    dt[, gene_id := trimws(gene_id)]
    
    # Aggregate reads by gene
    agg_dt <- dt[, .(Reads = sum(TPM)), by = gene_id]
    
    # Extract sample name safely based on generic logic
    sample_name <- sub("\\..*$", "", basename(file))
    setnames(agg_dt, "Reads", sample_name)
    return(agg_dt)
  }
  
  processed_list <- lapply(file_list, read_and_format)
  
  # Merge all data.tables
  merged_df <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), processed_list)
  
  # Convert to standard matrix
  mat <- as.matrix(merged_df[, -1, with = FALSE])
  rownames(mat) <- merged_df$gene_id
  mat[is.na(mat)] <- 0
  
  return(mat)
}

# ------------------------------------------------------------------------------
# 2. Translational Efficiency (TE) Calculation & Plotting
# ------------------------------------------------------------------------------
calculate_and_plot_te <- function(rna_mat, ribo_mat, target_gene, stage_vector) {
  
  # Validate Input
  if(!(target_gene %in% rownames(rna_mat)) | !(target_gene %in% rownames(ribo_mat))) {
    stop("Target gene not found in either RNA or Ribo matrices.")
  }
  
  # Extract vectors
  rna_expr <- as.numeric(rna_mat[target_gene, ])
  ribo_expr <- as.numeric(ribo_mat[target_gene, ])
  
  df <- data.frame(
    RNA = rna_expr,
    RPF = ribo_expr,
    Stage = factor(stage_vector, levels = unique(stage_vector))
  )
  
  # Calculate TE (add pseudocount to avoid division by zero)
  df$TE <- (df$RPF + 0.01) / (df$RNA + 0.01)
  
  # Scatter Plot
  p_scatter <- ggplot(df, aes(x = RNA, y = RPF, color = Stage)) +
    geom_point(size = 5) +
    scale_color_manual(values = met.brewer("Manet", length(unique(stage_vector)))) +
    labs(title = paste(target_gene, "Translation Dynamics"),
         x = "Normalized RNA Counts", y = "Normalized RPF Counts") +
    theme_science()
  
  # TE Trend line (Loess smoothing)
  p_trend <- ggplot(df, aes(x = Stage, y = TE, group = 1)) +
    geom_smooth(method = "loess", se = TRUE, color = "blue", linewidth = 1.2, alpha = 0.2) +
    geom_point(aes(color = Stage), size = 3) +
    scale_color_manual(values = met.brewer("Manet", length(unique(stage_vector)))) +
    labs(title = paste(target_gene, "Translational Efficiency (TE)"), x = "Developmental Stage", y = "TE (RPF/RNA)") +
    theme_science()
  
  return(list(Data = df, Scatter = p_scatter, Trend = p_trend))
}

# ------------------------------------------------------------------------------
# Example Execution
# ------------------------------------------------------------------------------
# rna_matrix <- merge_feature_files("path/to/RNA")
# ribo_matrix <- merge_feature_files("path/to/Ribo")
# stages <- c(rep('Zygote',2), rep('2C',2), rep('4C',2), rep('Morula',2), rep('BL',2))
# te_res <- calculate_and_plot_te(rna_matrix, ribo_matrix, "ENSMUSGXXXXXX", stages)
# print(te_res$Scatter)