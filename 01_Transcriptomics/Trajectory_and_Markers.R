# ==============================================================================
# Script: Trajectory_and_Markers.R
# Purpose: Cluster redefining, Monocle2 trajectory inference, and marker visualization.
# ==============================================================================
source("00_utils.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle)
  library(MetBrewer)
})

# ------------------------------------------------------------------------------
# 1. Trajectory Inference (Monocle 2)
# ------------------------------------------------------------------------------
run_monocle2_trajectory <- function(seurat_obj, group_by_col = "Celltype10", top_n_var = 1000) {
  
  # Ensure Idents are set correctly
  Idents(seurat_obj) <- group_by_col
  
  # Extract components for Monocle
  expr_matrix <- as(as.matrix(seurat_obj@assays$integrated@data), 'sparseMatrix')
  p_data <- new('AnnotatedDataFrame', data = seurat_obj@meta.data)
  f_data <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(expr_matrix), 
                                                        row.names = row.names(expr_matrix)))
  
  # Construct CellDataSet (assuming data is already log-normalized)
  cds <- newCellDataSet(expr_matrix, phenoData = p_data, featureData = f_data, expressionFamily = uninormal())
  
  # Feature selection
  var_genes <- VariableFeatures(seurat_obj)
  if(length(var_genes) < top_n_var) warning("Variable features are fewer than requested.")
  var_genes <- head(var_genes, top_n_var)
  
  cds <- setOrderingFilter(cds, var_genes)
  
  # Dimensionality reduction & Trajectory inference
  cds <- reduceDimension(cds, norm_method="none", reduction_method="DDRTree", 
                         max_components=6, scaling=TRUE, verbose=TRUE, pseudo_expr=0)
  
  # Order cells (reverse = TRUE logic can be adjusted based on known root state)
  cds <- orderCells(cds, reverse = TRUE)
  return(cds)
}

# ------------------------------------------------------------------------------
# 2. Visualization Functions
# ------------------------------------------------------------------------------
plot_marker_dots <- function(seurat_obj, genes, group_by, title = "Marker Expression") {
  
  # Filter genes that actually exist in the dataset
  valid_genes <- intersect(genes, rownames(seurat_obj))
  if(length(valid_genes) == 0) stop("None of the specified genes found in the Seurat object.")
  
  DefaultAssay(seurat_obj) <- 'RNA'
  
  p <- DotPlot(seurat_obj, features = valid_genes, group.by = group_by) +
    RotatedAxis() + coord_flip() +
    scale_color_gradientn(
      values = seq(0, 1, 0.2),
      colours = met.brewer("VanGogh2", n = 5, type = "discrete", direction = -1)
    ) +
    ggtitle(title) +
    theme_science()
  
  return(p)
}

# ------------------------------------------------------------------------------
# Example Execution
# ------------------------------------------------------------------------------
# 1. Monocle2 execution
# cds_obj <- run_monocle2_trajectory(embryo.integrated3, "Stage", 950)
# p_traj <- plot_cell_trajectory(cds_obj, color_by = "Stage", cell_size = 2) + 
#           scale_color_manual(values = opts$colors$stage_original) + theme_science()
# print(p_traj)

# 2. Marker visualization
# target_genes <- c('Pou5f1', 'Sox2', 'Nanog', 'Gata4', 'Gata6')
# print(plot_marker_dots(embryo.integrated3, target_genes, "Celltype10", "Pluripotency & Endoderm"))