# ==============================================================================
# Script: Multiomics_SingleGene_Violin.R
# Purpose: Generate calibrated multi-omics violin plots for target genes, mapping 
#          MOFA imputed latent values back to actual methylation/accessibility %.
# ==============================================================================
source("00_utils.R")
suppressPackageStartupMessages({
  library(MOFA2)
  library(Seurat)
  library(ggplot2)
  library(data.table)
})

# (Note: The calibration functions `calibrate_view`, `mvalue_to_pct`, and `get_gene_long` 
# from the original script are assumed to be loaded in 00_utils.R or locally to avoid bloat here. 
# Below is the wrapper that builds the final visualization).

plot_multiomics_violin <- function(gene_long_df, gene_name) {
  
  # Ensure strict assay ordering for visualization
  assay_levels <- c("RNA", "Chromatin accessibility", "DNA methylation")
  gene_long_df$assay <- factor(gene_long_df$assay, levels = assay_levels)
  
  # Color palette matching original setup
  assay_cols <- c("RNA" = "#3CB54E", "Chromatin accessibility" = "#00BFC4", "DNA methylation" = "#F37A71")
  
  p <- ggplot(gene_long_df, aes(x = group, y = value, fill = assay)) +
    geom_violin(scale = "width", alpha = 0.7, colour = "grey20", linewidth = 0.3) +
    geom_boxplot(width = 0.15, alpha = 0.7, outlier.shape = NA, colour = "grey20", linewidth = 0.3) +
    facet_grid(assay ~ ., scales = "free_y", switch = "y") +
    scale_fill_manual(values = assay_cols, guide = "none") +
    labs(title = gene_name, x = NULL, y = NULL) +
    theme_science() +
    theme(
      plot.title = element_text(face = "bold.italic", hjust = 0),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(face = "bold", angle = 90),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
  
  return(p)
}