# ==============================================================================
# Script: Seurat_DEG_and_Enrichment.R
# Purpose: Generalized Seurat FindAllMarkers/FindMarkers pipeline, Volcano plotting, 
#          and clusterProfiler GO/KEGG enrichment visualization.
# ==============================================================================
source("00_utils.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggrepel)
  library(viridis)
})

# ------------------------------------------------------------------------------
# 1. Standard Volcano Plot for Seurat Markers
# ------------------------------------------------------------------------------
plot_seurat_volcano <- function(marker_df, p_val_cut = 0.05, logfc_cut = 1, target_genes = NULL, title = "Volcano Plot") {
  
  # Categorize Significance
  marker_df$Sig <- "None"
  marker_df$Sig[marker_df$p_val < p_val_cut & marker_df$avg_log2FC >= logfc_cut] <- "Up"
  marker_df$Sig[marker_df$p_val < p_val_cut & marker_df$avg_log2FC <= -logfc_cut] <- "Down"
  marker_df$Sig <- factor(marker_df$Sig, levels = c("Up", "Down", "None"))
  
  # Handle labels
  marker_df$label <- ""
  if(!is.null(target_genes)) {
    hits <- intersect(target_genes, rownames(marker_df))
    marker_df[hits, "label"] <- hits
  }
  
  p <- ggplot(marker_df, aes(x = avg_log2FC, y = -log10(p_val), colour = Sig)) +
    geom_point(alpha = 0.6, size = 2.5) +
    scale_color_manual(values = c("Up" = "#ff4757", "Down" = "#546de5", "None" = "#d2dae2"), drop = FALSE) +
    geom_vline(xintercept = c(-logfc_cut, logfc_cut), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(p_val_cut), linetype = "dashed", color = "black") +
    labs(x = "log2(Fold Change)", y = "-log10 (P-value)", title = title) +
    theme_science() +
    theme(legend.position = "none")
  
  if(!is.null(target_genes) && length(hits) > 0) {
    p <- p + geom_label_repel(aes(label = label), size = 4, color = 'black', max.overlaps = Inf, fill = "white", alpha = 0.8)
  }
  
  return(p)
}

# ------------------------------------------------------------------------------
# 2. Batch GO/KEGG Enrichment
# ------------------------------------------------------------------------------
run_batch_enrichment <- function(marker_list, logfc_cut = 0.5, p_val_cut = 0.05, org_db = "org.Mm.eg.db", organism = "mmu") {
  
  res_list <- list(Up_GO = list(), Down_GO = list(), Up_KEGG = list(), Down_KEGG = list())
  
  for (name in names(marker_list)) {
    df <- marker_list[[name]]
    
    # Filter
    up_genes <- rownames(df[df$avg_log2FC >= logfc_cut & df$p_val < p_val_cut, ])
    down_genes <- rownames(df[df$avg_log2FC <= -logfc_cut & df$p_val < p_val_cut, ])
    
    # Translate IDs
    if(length(up_genes) > 0) {
      eg_up <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)$ENTREZID
      res_list$Up_GO[[name]] <- enrichGO(gene = eg_up, OrgDb = org_db, ont = "BP", pvalueCutoff = 0.05, readable = TRUE)
      res_list$Up_KEGG[[name]] <- enrichKEGG(gene = eg_up, organism = organism, pvalueCutoff = 0.05)
    }
    
    if(length(down_genes) > 0) {
      eg_down <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)$ENTREZID
      res_list$Down_GO[[name]] <- enrichGO(gene = eg_down, OrgDb = org_db, ont = "BP", pvalueCutoff = 0.05, readable = TRUE)
      res_list$Down_KEGG[[name]] <- enrichKEGG(gene = eg_down, organism = organism, pvalueCutoff = 0.05)
    }
  }
  return(res_list)
}

# ------------------------------------------------------------------------------
# 3. Generalized GO Bubble Plot
# ------------------------------------------------------------------------------
plot_go_bubble <- function(go_df, group_col = "Cluster", term_col = "Description", pval_col = "pvalue", count_col = "Count") {
  
  # Ensure factors are set if user wants specific ordering
  if(is.character(go_df[[term_col]])) {
    go_df[[term_col]] <- factor(go_df[[term_col]], levels = unique(go_df[[term_col]]))
  }
  
  p <- ggplot(go_df) +
    geom_point(aes(x = .data[[group_col]], y = .data[[term_col]], color = -log10(.data[[pval_col]]), size = .data[[count_col]])) +
    scale_color_gradientn(colors = rev(magma(256))) +
    scale_size_continuous(range = c(2, 8)) +
    theme_science() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Group", y = "Enriched Term", color = "-log10(P-val)", size = "Gene Count")
  
  return(p)
}

# ------------------------------------------------------------------------------
# Example Execution
# ------------------------------------------------------------------------------
# Idents(embryo.integrated3) <- 'Celltype11'
# markers <- FindMarkers(embryo.integrated3, ident.1="ClusterA", ident.2="ClusterB", test.use="MAST", logfc.threshold=1)
# target_list <- c("Gene1", "Gene2", "Gene3")
# print(plot_seurat_volcano(markers, target_genes = target_list, title = "ClusterA vs ClusterB"))
#
# go_res <- run_batch_enrichment(list(Comp1 = markers))