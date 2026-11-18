# ==============================================================================
# Script: MultiOmics_Correlation.R
# Purpose: Calculate and visualize correlation between chromatin accessibility (acc), 
#          DNA methylation (met), and RNA expression.
# ==============================================================================
source("00_utils.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(colorspace)
})

# ------------------------------------------------------------------------------
# 1. 2D Density Plot for Sc-Multiomics (e.g. Acc vs RNA)
# ------------------------------------------------------------------------------
plot_multiomics_density <- function(df, x_col, y_col, x_label = "Accessibility Rate", y_label = "RNA Expression") {
  
  # Ensure columns are numeric
  df[[x_col]] <- as.numeric(df[[x_col]])
  df[[y_col]] <- as.numeric(df[[y_col]])
  
  # Calculate Correlation
  cor_res <- cor.test(df[[x_col]], df[[y_col]], method = "spearman")
  r_val <- round(cor_res$estimate, 2)
  p_val <- signif(cor_res$p.value, 2)
  
  # Plot
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
    scale_fill_continuous_sequential(palette = "Rocket") +
    geom_point(color = "black", size = 0.5, alpha = 0.3) +
    labs(x = x_label, y = y_label, 
         title = "Multi-omics Correlation",
         subtitle = paste0("Spearman R = ", r_val, " (p = ", p_val, ")")) +
    theme_science() +
    theme(legend.position = "right")
  
  return(p)
}

# ------------------------------------------------------------------------------
# 2. Smooth Scatter Plot for Replicate Correlation (Red Palette)
# ------------------------------------------------------------------------------
plot_replicate_smooth_scatter <- function(raw_counts_file, title = "Replicate Correlation") {
  
  df <- read.delim(raw_counts_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
  
  # Assuming deepTools multiBigwigSummary outputs counts in columns 4 and 5
  x_raw <- as.numeric(df[, 4])
  y_raw <- as.numeric(df[, 5])
  
  valid <- is.finite(x_raw) & is.finite(y_raw)
  keep <- (x_raw[valid] > 0 | y_raw[valid] > 0)
  x_clean <- x_raw[valid][keep]
  y_clean <- y_raw[valid][keep]
  
  r_val <- ifelse(length(x_clean) > 2, cor(x_clean, y_clean, method = "pearson"), NA)
  
  x_plot <- log2(x_clean + 1)
  y_plot <- log2(y_clean + 1)
  
  reds <- colorRampPalette(c("#fff5f0", "#fcbba1", "#fb6a4a", "#cb181d", "#67000d"))
  
  # Open graphics device if saving is needed, here we just plot to active device
  par(mar = c(4.5, 4.5, 3, 1.5))
  smoothScatter(x_plot, y_plot, colramp = reds, nbin = 128, nrpoints = 150, pch = 20, cex = 0.4, 
                col = "black", xlab = "Rep1 log2 (RPKM+1)", ylab = "Rep2 log2 (RPKM+1)", 
                main = title, bty = "l", cex.axis = 1.2, cex.lab = 1.3, cex.main = 1.6)
  
  if (!is.na(r_val)) {
    legend("topleft", bty = "n", legend = paste0("R=", sprintf("%.2f", r_val)), cex = 1.5, text.font = 2)
  }
}