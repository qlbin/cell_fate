# ==============================================================================
# Script: hdWGCNA_Analysis.R
# Purpose: Generalized High-Dimensional WGCNA (hdWGCNA) pipeline for scRNA-seq.
# ==============================================================================
source("00_utils.R")
suppressPackageStartupMessages({
  library(WGCNA)
  library(hdWGCNA)
  library(UCell)
  library(MetBrewer)
  library(GeneOverlap)
  library(enrichR)
  library(viridis)
  library(ggrepel)
})

# ------------------------------------------------------------------------------
# 1. Initialize and Construct hdWGCNA Network
# ------------------------------------------------------------------------------
run_hdwgcna_pipeline <- function(seurat_obj, group_by_col, selected_groups = NULL, 
                                 wgcna_name = "scWGCNA", fraction = 0.05, soft_power = 8, min_module_size = 100) {
  
  DefaultAssay(seurat_obj) <- 'RNA'
  
  # Setup
  seurat_obj <- SetupForWGCNA(seurat_obj, gene_select = "fraction", fraction = fraction, wgcna_name = wgcna_name)
  
  # Construct Metacells
  seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj, group.by = group_by_col, reduction = 'umap',
    k = 5, max_shared = 15, min_cells = 5, ident.group = group_by_col
  )
  seurat_obj <- NormalizeMetacells(seurat_obj)
  
  # Set DatExpr
  if (is.null(selected_groups)) {
    selected_groups <- unique(seurat_obj@meta.data[[group_by_col]])
  }
  seurat_obj <- SetDatExpr(seurat_obj, group_name = selected_groups, group.by = group_by_col, assay = 'RNA', slot = 'data')
  
  # Network Construction (Assumes TestSoftPowers has already been evaluated by user; forcing soft_power)
  seurat_obj <- ConstructNetwork(seurat_obj, soft_power = soft_power, overwrite_tom = TRUE, 
                                 group.by = group_by_col, minModuleSize = min_module_size)
  
  # Scale data and Compute MEs
  seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
  seurat_obj <- ModuleEigengenes(seurat_obj, group.by.vars = group_by_col)
  seurat_obj <- ModuleConnectivity(seurat_obj, harmonized = FALSE)
  
  # Rename modules and adjust colors systematically
  seurat_obj <- ResetModuleNames(seurat_obj, new_name = "M")
  n_mods <- nrow(GetModules(seurat_obj) %>% dplyr::select(c(module, color)) %>% distinct()) - 1
  new_colors <- paste0(met.brewer("Signac", n = n_mods))
  seurat_obj <- ResetModuleColors(seurat_obj, new_colors)
  
  return(seurat_obj)
}

# ------------------------------------------------------------------------------
# 2. Extract Hubs & Score Cells
# ------------------------------------------------------------------------------
score_wgcna_modules <- function(seurat_obj, group_by_col, n_hubs = 100) {
  
  hub_df <- GetHubGenes(seurat_obj, n_hubs = n_hubs)
  # Optional: Convert to ENTREZID if needed
  # hub_df$ENTREZID <- bitr(hub_df$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb='org.Mm.eg.db', drop=FALSE)$ENTREZID
  
  # Score cells based on module hubs
  seurat_obj <- ModuleExprScore(seurat_obj, n_genes = n_hubs, method = 'Seurat')
  
  # Append MEs to metadata for Seurat plotting
  MEs <- GetMEs(seurat_obj, harmonized = FALSE)
  mods <- colnames(MEs)[colnames(MEs) != 'grey']
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
  
  p_dot <- DotPlot(seurat_obj, features = mods, group.by = group_by_col) +
    scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
    RotatedAxis() + theme_science()
  
  return(list(Object = seurat_obj, Hubs = hub_df, DotPlot = p_dot))
}

# ------------------------------------------------------------------------------
# 3. Module vs Marker Overlap
# ------------------------------------------------------------------------------
plot_module_marker_overlap <- function(seurat_obj, marker_df, fc_cutoff = 1) {
  
  overlap_df <- OverlapModulesDEGs(seurat_obj, deg_df = marker_df, fc_cutoff = fc_cutoff)
  
  p_overlap <- OverlapDotPlot(overlap_df, plot_var = 'odds_ratio') +
    ggtitle('Overlap of WGCNA Modules & Cell-Type Markers') + theme_science()
  
  return(list(Data = overlap_df, Plot = p_overlap))
}

# ------------------------------------------------------------------------------
# Example Execution
# ------------------------------------------------------------------------------
# wgcna_res <- run_hdwgcna_pipeline(embryo.integrated3, "Celltype10")
# scored_res <- score_wgcna_modules(wgcna_res, "Celltype10")
# print(scored_res$DotPlot)
# 
# markers <- Seurat::FindAllMarkers(wgcna_res, only.pos=TRUE, logfc.threshold=1)
# overlap_res <- plot_module_marker_overlap(wgcna_res, markers)
# print(overlap_res$Plot)