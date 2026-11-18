# ==============================================================================
# Script: External_Data_Integration.R
# Purpose: Process external single-cell RNA-seq datasets and perform advanced 
#          expression visualization (e.g., pseudo-log plots for early development).
# ==============================================================================
source("00_utils.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(scales)
  library(tidyr)
  library(MetBrewer)
})

# ------------------------------------------------------------------------------
# 1. Load and Merge External Matrices
# ------------------------------------------------------------------------------
build_external_seurat <- function(matrix_dir, file_pattern = "\\.txt\\.gz$") {
  
  file_list <- list.files(path = matrix_dir, pattern = file_pattern, full.names = TRUE)
  if(length(file_list) == 0) stop("No matrix files found.")
  
  read_and_format <- function(file) {
    dt <- fread(file, check.names = TRUE)
    
    # Handle variations in gene symbol column names
    if ("X.Gene_symbol" %in% colnames(dt)) setnames(dt, "X.Gene_symbol", "Gene_symbol")
    if (!"Gene_symbol" %in% colnames(dt)) setnames(dt, 1, "Gene_symbol")
    
    dt[, Gene_symbol := trimws(Gene_symbol)]
    
    # Aggregate duplicate symbols
    if("reads" %in% colnames(dt)) {
      unique_dt <- dt[, .(Reads = sum(reads)), by = Gene_symbol]
    } else {
      # Assuming 2nd column contains counts if not named 'reads'
      val_col <- colnames(dt)[2]
      unique_dt <- dt[, .(Reads = sum(get(val_col))), by = Gene_symbol]
    }
    
    # Extract sample name from file (Adjust the split logic as per your file naming convention)
    file_name <- basename(file)
    parts <- strsplit(file_name, "_")[[1]]
    if(length(parts) >= 3) {
      new_name <- paste(parts[2], parts[3], sep = "_")
    } else {
      new_name <- sub(file_pattern, "", file_name)
    }
    
    setnames(unique_dt, "Reads", new_name)
    return(unique_dt)
  }
  
  processed_list <- lapply(file_list, read_and_format)
  final_df <- Reduce(function(x, y) merge(x, y, by = "Gene_symbol", all = TRUE), processed_list)
  
  mat <- as.matrix(final_df[, -1, with = FALSE])
  rownames(mat) <- final_df$Gene_symbol
  mat[is.na(mat)] <- 0
  
  # Standard Seurat Processing
  obj <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 2000)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 5000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = min(50, ncol(obj)-1), verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:min(30, ncol(obj)-1))
  obj <- FindClusters(obj, resolution = 1.0)
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:min(30, ncol(obj)-1))
  
  return(obj)
}

# ------------------------------------------------------------------------------
# 2. Advanced Expression Visualization (Pseudo-log scale)
# ------------------------------------------------------------------------------
#' Plot gene expression dynamics using pseudo-log scale (optimal for low-expression genes)
plot_pseudolog_dynamics <- function(seurat_obj, genes, group_by_col, cell_order = NULL, sigma = 0.01) {
  
  # Fetch Data
  expr_data <- FetchData(seurat_obj, vars = c(genes, group_by_col))
  colnames(expr_data)[ncol(expr_data)] <- "Group"
  expr_data$CellID <- rownames(expr_data)
  
  if(!is.null(cell_order)) {
    expr_data$Group <- factor(expr_data$Group, levels = cell_order)
    expr_data <- expr_data[!is.na(expr_data$Group), ]
  }
  
  gene_long <- expr_data %>%
    pivot_longer(cols = all_of(genes), names_to = "Gene", values_to = "Expression")
  
  # Apply Pseudo-log transformation
  pseudo_trans <- pseudo_log_trans(sigma = sigma)
  gene_long$Expr_pseudolog <- pseudo_trans$transform(gene_long$Expression)
  
  # Calculate means in the transformed space
  mean_expr <- gene_long %>%
    group_by(Gene, Group) %>%
    summarise(Mean_pseudolog = mean(Expr_pseudolog, na.rm = TRUE), .groups = "drop")
  
  # Plotting
  p <- ggplot(gene_long, aes(x = Group, y = Expr_pseudolog)) +
    geom_dotplot(binaxis = "y", stackdir = "center", fill = "lightgreen", dotsize = 0.1, shape = 21, color = NA) +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange", shape = 21, fill = "white", size = 0.3) +
    geom_line(data = mean_expr, aes(x = Group, y = Mean_pseudolog, group = 1), color = "darkgreen", linewidth = 1) +
    facet_wrap(~ Gene, scales = "free_y") +
    labs(x = "Cell Type / Group", y = paste0("Expression Level (pseudo-log, \u03C3=", sigma, ")")) +
    theme_science() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# ------------------------------------------------------------------------------
# Example Execution
# ------------------------------------------------------------------------------
# ext_obj <- build_external_seurat("path/to/matrix/folder")
# target_genes <- c("GeneA", "GeneB")
# stages <- c("Zygote", "2C", "4C", "8C", "16C", "32C", "64C", "BL")
# print(plot_pseudolog_dynamics(ext_obj, target_genes, "orig.ident", stages))