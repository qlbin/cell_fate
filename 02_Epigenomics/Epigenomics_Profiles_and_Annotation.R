# ==============================================================================
# Script: Epigenomics_Profiles_and_Annotation.R
# Purpose: Advanced plotting of deepTools profiles (baseline-subtracted) 
#          and ChIPseeker genomic annotation.
# ==============================================================================
source("00_utils.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(ChIPseeker)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(TxDb.Mmusculus.UCSC.mm39.knownGene)
})

# ------------------------------------------------------------------------------
# 1. Advanced deepTools Profile Plotting (Baseline Subtracted)
# ------------------------------------------------------------------------------
plot_deeptools_profile <- function(file_path, sample_levels, colors, exclude_samples = NULL, 
                                   baseline_prob = 0.05, title = "Epigenomic Signal Profile") {
  
  # Read deepTools raw matrix (header=FALSE)
  raw_df <- read.delim(file_path, header = FALSE, stringsAsFactors = FALSE)
  
  # Extract relevant rows based on provided sample names
  target_rows <- raw_df[raw_df[, 1] %in% sample_levels, ]
  rownames(target_rows) <- target_rows[, 1]
  
  # Transpose and clean
  df_clean <- as.data.frame(t(target_rows[, -1]))
  df_clean <- df_clean[, intersect(sample_levels, colnames(df_clean))]
  suppressWarnings(df_clean[] <- lapply(df_clean, as.numeric))
  df_clean <- df_clean[complete.cases(df_clean), ]
  
  # Defensive truncation (e.g., if multiple groups were stacked by deepTools)
  if (nrow(df_clean) > 400) { df_clean <- df_clean[1:(nrow(df_clean) / 2), ] }
  
  # Generate X-axis (assuming standard -3kb to 3kb)
  df_clean$Distance <- seq(-3000, 3000, length.out = nrow(df_clean))
  
  df_long <- pivot_longer(df_clean, cols = -Distance, names_to = "Sample", values_to = "Signal")
  df_long$Sample <- factor(df_long$Sample, levels = sample_levels)
  
  # Baseline Subtraction
  df_long <- df_long %>%
    group_by(Sample) %>%
    mutate(
      baseline = quantile(Signal, probs = baseline_prob, na.rm = TRUE),
      Signal_norm = pmax(Signal - baseline, 0)
    ) %>%
    ungroup()
  
  # Optional Exclusion
  if (!is.null(exclude_samples)) {
    df_long <- subset(df_long, !(Sample %in% exclude_samples))
  }
  
  p <- ggplot(df_long, aes(x = Distance, y = Signal_norm, color = Sample)) +
    geom_line(linewidth = 0.8, alpha = 0.9) +  
    scale_color_manual(values = colors) +
    labs(title = title, x = "Distance from center (bp)", y = "Signal (Baseline Subtracted)") +
    scale_x_continuous(breaks = c(-3000, 0, 3000), labels = c("-3kb", "center", "3kb"), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_science() + 
    theme(legend.position = c(0.85, 0.75), legend.background = element_rect(fill = alpha("white", 0.7), color = NA))
  
  return(p)
}

# ------------------------------------------------------------------------------
# 2. ChIPseeker Annotation & Enrichment wrapper
# ------------------------------------------------------------------------------
annotate_and_enrich_peaks <- function(peak_list, txdb = TxDb.Mmusculus.UCSC.mm39.knownGene, orgDb = "org.Mm.eg.db") {
  
  res_list <- list()
  
  # Annotate Peaks
  anno_list <- lapply(peak_list, annotatePeak, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = orgDb, verbose = FALSE)
  res_list$Annotations <- anno_list
  
  # Extract Gene Symbols
  genes_symbol <- lapply(anno_list, function(x) {
    gene_ids <- as.data.frame(x)$geneId
    bitr(gene_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = orgDb, drop = FALSE)
  })
  res_list$Symbols <- genes_symbol
  
  # Optional: GO Enrichment if names are provided
  if (!is.null(names(peak_list))) {
    entrez_list <- lapply(anno_list, function(x) as.data.frame(x)$geneId)
    compGO <- compareCluster(geneCluster = entrez_list, OrgDb = orgDb, fun = "enrichGO", pvalueCutoff = 0.05)
    res_list$GO_Plot <- dotplot(compGO, showCategory = 15, title = "Peak GO Enrichment") + theme_science()
  }
  
  return(res_list)
}

# ------------------------------------------------------------------------------
# Example Execution
# ------------------------------------------------------------------------------
# my_colors <- c("Ctrl" = "#7F7F7F", "mt" = "#00617F", "NLS" = "#A80735")
# p_profile <- plot_deeptools_profile("profile_data.txt", c("Ctrl", "mt", "NLS"), my_colors)
# print(p_profile)