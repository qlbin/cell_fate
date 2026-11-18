# ==============================================================================
# Script: SCENIC_Regulon_Analysis.R
# Purpose: Run SCENIC pipeline, evaluate Regulon Specificity Score (RSS), 
#          and perform functional enrichment on Regulon target genes.
# ==============================================================================
source("00_utils.R")
suppressPackageStartupMessages({
  library(SCENIC)
  library(Seurat)
  library(AUCell)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(UpSetR)
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 1. SCENIC Pipeline Wrapper (Assume scenicOptions is pre-configured)
# ------------------------------------------------------------------------------
run_scenic_wrapper <- function(seurat_obj, scenicOptions) {
  
  DefaultAssay(seurat_obj) <- 'RNA'
  exprMat <- as.matrix(GetAssayData(seurat_obj, slot = "counts"))
  
  # Filter genes
  motifRankings <- importRankings(getDatabases(scenicOptions)[[1]])
  genesKept <- intersect(rownames(exprMat), colnames(getRanking(motifRankings)))
  exprMat_filtered <- exprMat[genesKept, ]
  
  # Correlation and GENIE3
  runCorrelation(exprMat_filtered, scenicOptions)
  runGenie3(log2(exprMat_filtered + 1), scenicOptions)
  
  # Network to Regulons
  scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
  scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
  scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, log2(exprMat + 1))
  scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
  
  return(scenicOptions)
}

# ------------------------------------------------------------------------------
# 2. Regulon Specificity Score (RSS) & Visualization
# ------------------------------------------------------------------------------
plot_top_rss <- function(scenicOptions, seurat_obj, group_by_col, top_n = 5) {
  
  regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
  cellInfo <- seurat_obj@meta.data
  
  # Calculate RSS
  rss <- calcRSS(AUC = getAUC(regulonAUC), cellAnnotation = cellInfo[[group_by_col]])
  rssPlot <- plotRSS(rss)
  
  rss_data <- rssPlot$df
  colnames(rss_data) <- c("Regulon", "Group", "RSS", "Z", "Number")
  
  # Extract top N regulons per group
  top_rss <- rss_data %>% group_by(Group) %>% top_n(n = top_n, wt = RSS) %>% as.data.frame()
  
  # Clean TF names for expression plotting
  top_rss$TF_Gene <- stringr::str_extract(top_rss$Regulon, "[^ _]+")
  
  # Fetch Expression Data
  p_dot <- DotPlot(seurat_obj, features = unique(top_rss$TF_Gene), group.by = group_by_col)
  expr_data <- p_dot$data
  
  # Merge RSS and Expression
  top_rss$combined <- paste(top_rss$Group, top_rss$TF_Gene, sep = "_")
  expr_data$combined <- paste(expr_data$id, expr_data$features.plot, sep = "_")
  
  merged_df <- merge(top_rss, expr_data, by = "combined")
  merged_df <- na.omit(merged_df)
  
  p_final <- ggplot(merged_df) +
    geom_point(aes(x = Group, y = reorder(features.plot, Group), color = avg.exp.scaled, size = RSS)) +
    scale_color_gradientn(values = seq(0, 1, 0.2), colours = c("lightgrey", "#4E68A3", "#DDB745", "#E67336", "#DA2F2F")) +
    theme_science() + RotatedAxis() +
    labs(x = "Cell Group", y = "Top Regulon TFs", title = "Regulon Specificity & TF Expression")
  
  return(list(RSS_Data = rss_data, Top_RSS = top_rss, Plot = p_final))
}

# ------------------------------------------------------------------------------
# 3. Regulon Target Gene Enrichment & UpSet Plot
# ------------------------------------------------------------------------------
enrich_regulon_targets <- function(scenicOptions, target_tfs) {
  
  regulons <- loadInt(scenicOptions, "regulons")
  matched_regulons <- regulons[intersect(names(regulons), target_tfs)]
  
  if(length(matched_regulons) == 0) stop("None of the target TFs found in Regulons.")
  
  geneList_entrez <- lapply(matched_regulons, function(genes) {
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
  })
  
  go_results <- compareCluster(geneCluster = geneList_entrez, fun = "enrichGO", OrgDb = org.Mm.eg.db, ont = "BP", pvalueCutoff = 0.05)
  
  # Prepare UpSet data
  go_df <- as.data.frame(go_results)
  if(nrow(go_df) > 0) {
    tf_go_list <- go_df %>% filter(p.adjust < 0.05) %>% group_by(Cluster) %>% summarise(GO_terms = list(ID)) %>% tibble::deframe()
    all_terms <- unique(unlist(tf_go_list))
    
    binary_mat <- sapply(tf_go_list, function(go_ids) all_terms %in% go_ids) * 1
    rownames(binary_mat) <- all_terms
    binary_df <- as.data.frame(binary_mat)
    
    upset_plot <- upset(binary_df, nsets = ncol(binary_df), nintersects = 20, order.by = "freq", sets.bar.color = "#2c7fb8")
    return(list(GO_Results = go_df, UpSet = upset_plot))
  } else {
    warning("No significant GO terms found for these Regulons.")
    return(NULL)
  }
}