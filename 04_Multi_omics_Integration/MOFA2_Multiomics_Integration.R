# ==============================================================================
# Script: MOFA2_Multiomics_Integration.R
# Purpose: Prepare multi-omics matrices, train MOFA2 model, and plot variance.
# ==============================================================================
source("00_utils.R")
suppressPackageStartupMessages({
  library(MOFA2)
  library(data.table)
  library(dplyr)
  library(reticulate)
})

# Use reticulate to bind to your specific python/conda environment before running
# use_python("/path/to/python", required=TRUE)

# ------------------------------------------------------------------------------
# 1. MOFA Model Training Wrapper
# ------------------------------------------------------------------------------
train_mofa_multiomics <- function(matrix_list, meta_df, group_col, outfile = "mofa_model.hdf5") {
  
  # Create Object
  mofa_obj <- create_mofa_from_matrix(matrix_list, groups = meta_df[[group_col]])
  
  # Set Options
  data_opts <- get_default_data_options(mofa_obj)
  model_opts <- get_default_model_options(mofa_obj)
  train_opts <- get_default_training_options(mofa_obj)
  
  mofa_obj <- prepare_mofa(mofa_obj, data_options = data_opts, model_options = model_opts, training_options = train_opts)
  
  # Train
  mofa_trained <- run_mofa(mofa_obj, outfile)
  
  # Append full metadata
  samples_metadata(mofa_trained) <- meta_df
  return(mofa_trained)
}

# ------------------------------------------------------------------------------
# 2. Standard MOFA Visualizations
# ------------------------------------------------------------------------------
plot_mofa_diagnostics <- function(mofa_model, group_col) {
  p_var <- plot_variance_explained(mofa_model, x = "view", y = "factor", split_by = "group") + RotatedAxis() + theme_science()
  
  mofa_model <- run_umap(mofa_model, factors = 1:10, n_neighbors = 15, min_dist = 0.5)
  p_umap <- plot_dimred(mofa_model, method = "UMAP", color_by = group_col, dot_size = 2) + theme_science()
  
  return(list(Variance = p_var, UMAP = p_umap))
}

export_top_weights <- function(mofa_model, views, factors = 1:5, top_n = 200) {
  weight_list <- list()
  for (v in views) {
    for (f in factors) {
      w <- get_weights(mofa_model, views = v, factor = f)[[1]][,1]
      weight_list[[paste(v, "Factor", f, sep="_")]] <- names(tail(sort(abs(w)), n = top_n))
    }
  }
  return(weight_list)
}