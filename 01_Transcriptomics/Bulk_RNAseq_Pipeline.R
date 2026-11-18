# ==============================================================================
# Script: Bulk_RNAseq_Pipeline.R
# Purpose: Generalized Bulk RNA-seq processing (TPM, Batch Correction, DESeq2, GO/GSEA)
# ==============================================================================
source("00_utils.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggrepel)
})

# ------------------------------------------------------------------------------
# 1. Normalization & Batch Correction
# ------------------------------------------------------------------------------
normalize_and_correct <- function(count_matrix, gene_lengths, batch_info, conditions) {
  
  # TPM Calculation
  gene_length_kb <- gene_lengths / 1000
  rpk <- count_matrix / gene_length_kb
  rpk[is.na(rpk) | !is.finite(rpk)] <- 0
  
  scaling_factor <- colSums(rpk)
  tpm <- t(t(rpk) / scaling_factor) * 1e6
  
  # Log2 transformation
  log_tpm <- log2(tpm + 1)
  
  # Batch correction preserving biological conditions
  design_mat <- model.matrix(~ 0 + conditions)
  corrected_tpm <- removeBatchEffect(x = log_tpm, batch = batch_info, design = design_mat)
  
  return(list(TPM = tpm, LogTPM = log_tpm, Corrected_LogTPM = corrected_tpm))
}

# ------------------------------------------------------------------------------
# 2. DESeq2 Differential Expression
# ------------------------------------------------------------------------------
run_deseq2_analysis <- function(raw_counts, meta_df, condition_col, ctrl_level, trt_level, min_reads = 10, min_samples = 3) {
  
  # Pre-filtering
  keep <- rowSums(raw_counts >= min_reads) >= min_samples
  filtered_counts <- round(raw_counts[keep, ])
  storage.mode(filtered_counts) <- "integer"
  
  # Make sure condition is factored properly
  meta_df[[condition_col]] <- factor(meta_df[[condition_col]], levels = c(ctrl_level, trt_level))
  
  # Run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = filtered_counts, colData = meta_df, design = as.formula(paste0("~", condition_col)))
  dds <- DESeq(dds, quiet = TRUE)
  
  res <- results(dds, contrast = c(condition_col, trt_level, ctrl_level))
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  
  # Return complete sorted results
  return(res_df[order(res_df$pvalue), ])
}

# ------------------------------------------------------------------------------
# 3. Enrichment (GO / GSEA) & Visualization
# ------------------------------------------------------------------------------
run_go_enrichment <- function(genes, org_db = org.Mm.eg.db, keyType = "SYMBOL") {
  ego <- enrichGO(
    gene = genes, OrgDb = org_db, keyType = keyType, ont = "BP",
    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE
  )
  return(ego)
}

plot_generalized_volcano <- function(res_df, title, lfc_cut = 1, p_cut = 0.05, target_genes = NULL) {
  
  res_df$Significance <- "NS"
  res_df$Significance[which(res_df$pvalue < p_cut & res_df$log2FoldChange > lfc_cut)] <- "Up"
  res_df$Significance[which(res_df$pvalue < p_cut & res_df$log2FoldChange < -lfc_cut)] <- "Down"
  res_df$Significance <- factor(res_df$Significance, levels = c("Up", "Down", "NS"))
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = Significance), alpha = 0.8, size = 1.5) +
    scale_color_manual(values = c("Up" = "#e41a1c", "Down" = "#377eb8", "NS" = "grey80"), drop = FALSE) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(p_cut), linetype = "dashed", color = "black") +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10(P-value)") +
    theme_science()
  
  if(!is.null(target_genes)) {
    hits <- res_df[res_df$gene_id %in% target_genes, ]
    if(nrow(hits) > 0) {
      p <- p + geom_text_repel(data = hits, aes(label = gene_id), size = 5, 
                               fontface = "bold.italic", box.padding = 0.5)
    }
  }
  return(p)
}