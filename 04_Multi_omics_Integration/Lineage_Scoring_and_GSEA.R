# ==============================================================================
# Script: Lineage_Scoring_and_GSEA.R
# Purpose: Compute generic gene set scores (log of summed CPM) and perform 
#          Fast Gene Set Enrichment Analysis (fgsea).
# ==============================================================================
source("00_utils.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(fgsea)
  library(org.Mm.eg.db)
  library(dplyr)
  library(ggplot2)
  library(ggsignif)
})

# ------------------------------------------------------------------------------
# 1. Gene Set Scoring (log-sum of CPM)
# ------------------------------------------------------------------------------
compute_geneset_scores <- function(seurat_obj, geneset_list, group_by_col) {
  
  # Extract count matrix and compute CPM
  counts_mat <- as.matrix(GetAssayData(seurat_obj, slot = "counts"))
  cpm_mat <- edgeR::cpm(counts_mat)
  
  score_df <- data.frame(Sample = colnames(seurat_obj), 
                         Group = seurat_obj@meta.data[[group_by_col]])
  
  for (set_name in names(geneset_list)) {
    valid_genes <- intersect(geneset_list[[set_name]], rownames(cpm_mat))
    if (length(valid_genes) > 0) {
      # Log of summed CPMs
      score_df[[set_name]] <- log(colSums(cpm_mat[valid_genes, , drop = FALSE]) + 1)
    } else {
      warning(paste("No valid genes found for", set_name))
    }
  }
  return(score_df)
}

plot_score_boxplot <- function(score_df, score_col, comparisons_list, title = "Gene Set Score") {
  ggplot(score_df, aes(x = Group, y = .data[[score_col]], fill = Group)) + 
    geom_boxplot() +
    geom_signif(comparisons = comparisons_list, step_increase = 0.1, test = wilcox.test, map_signif_level = TRUE) +
    theme_science() + RotatedAxis() +
    labs(title = title, x = "Group", y = "Log(Summed CPM + 1)") +
    theme(legend.position = "none")
}

# ------------------------------------------------------------------------------
# 2. Generalized fgsea Pipeline
# ------------------------------------------------------------------------------
run_fgsea_pipeline <- function(seurat_obj, ident.1, ident.2, gmt