# ==============================================================================
# Script: Motif_Enrichment_and_Integration.R
# Purpose: Parse MEME (AME) / JASPAR motif enrichment from differential ATAC peaks 
#          and integrate with corresponding TF expression from scRNA-seq.
# ==============================================================================
source("00_utils.R")
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(Seurat)
  library(monaLisa)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Mmusculus.UCSC.mm39)
})

# ------------------------------------------------------------------------------
# 1. JASPAR Motif Enrichment (monaLisa) for DiffBind Peaks
# ------------------------------------------------------------------------------
run_jaspar_enrichment <- function(granges_peaks, bg_method = "random", top_n = 10) {
  
  # Fetch Vertebrate Motifs from JASPAR
  opts <- list(species = "Mus musculus", all_versions = FALSE)
  pwms <- getMatrixSet(JASPAR2020, opts)
  
  # Extract sequences
  seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm39, granges_peaks)
  
  # Run Enrichment
  se_enrich <- calcBinnedMotifEnrR(seqs = seqs, pwms = pwms, background = bg_method)
  
  # Plot
  p <- plotMotifEnrHeatmap(se_enrich, nMotifs = top_n)
  return(list(Enrichment = se_enrich, Heatmap = p))
}

# ------------------------------------------------------------------------------
# 2. Integrate MEME (AME) Output with scRNA-seq TF Expression
# ------------------------------------------------------------------------------
plot_ame_rna_integration <- function(ame_dir, seurat_obj, group_by_col, target_tfs = NULL, 
                                     pval_cutoff = 0.05, top_n_per_group = 10) {
  
  # Parse AME TSV files from subdirectories
  folders <- list.dirs(ame_dir, full.names = TRUE, recursive = FALSE)
  folders <- folders[basename(folders) != "fasta"]
  
  read_ame <- function(f) {
    lines <- readLines(f)
    lines <- lines[!grepl("^#", lines) & nzchar(lines)]
    if (length(lines) <= 1) return(NULL)
    fread(text = paste(lines, collapse = "\n"), sep = "\t")
  }
  
  ame_list <- lapply(folders, function(folder) {
    f <- file.path(folder, "ame.tsv")
    if (file.exists(f)) {
      dt <- read_ame(f)
      if (!is.null(dt)) {
        dt[, group := basename(folder)]
        return(dt)
      }
    }
    return(NULL)
  })
  
  ame <- rbindlist(ame_list, fill = TRUE)
  if(nrow(ame) == 0) stop("No valid AME results found.")
  
  # Clean Motif Names (e.g., POU5F1::SOX2 -> Pou5f1)
  ame <- ame[, .(motif = unlist(strsplit(motif_alt_ID, "::"))), by = setdiff(names(ame), "motif_alt_ID")]
  ame[, motif := sub("\\(.*\\)$", "", motif)]
  ame[, gene := paste0(toupper(substr(motif, 1, 1)), tolower(substring(motif, 2)))]
  ame[, log10p := -log10(as.numeric(`p-value`) + 1e-300)]
  
  # Filter Motifs
  if (!is.null(target_tfs)) {
    sel <- paste0(toupper(substr(target_tfs, 1, 1)), tolower(substring(target_tfs, 2)))
    plot_dt <- ame[gene %in% sel]
  } else {
    plot_dt <- ame[order(as.numeric(`p-value`))][, head(.SD, top_n_per_group), by = group]
  }
  
  # Deduplicate same TF in same group
  plot_dt <- plot_dt[order(as.numeric(`p-value`))][, .SD[1], by = .(group, gene)]
  
  # Fetch RNA Expression
  valid_genes <- intersect(unique(plot_dt$gene), rownames(seurat_obj))
  if(length(valid_genes) == 0) stop("None of the enriched TFs found in Seurat object.")
  
  p_dot <- DotPlot(seurat_obj, features = valid_genes, group.by = group_by_col)
  expr_data <- as.data.table(p_dot$data)
  
  # Standardize integration keys
  # (Assuming your AME group names contain stage info that maps to Seurat Idents)
  # Here we use a generic mapping approach; adjust substring logic based on exact AME folder naming
  plot_dt[, stage_up := gsub(".*_up_in_", "", group)] 
  expr_data[, `:=`(gene = as.character(features.plot), stage_up = as.character(id))]
  
  final_dt <- merge(plot_dt, expr_data, by = c("gene", "stage_up"), all.x = TRUE)
  final_dt[is.na(avg.exp.scaled), avg.exp.scaled := 0]
  
  # Cap sizes for aesthetics
  cap <- 50 
  final_dt[, size_val := pmin(log10p, cap)]
  
  p <- ggplot(final_dt, aes(x = stage_up, y = reorder(gene, log10p))) +
    geom_point(aes(size = size_val, fill = avg.exp.scaled), shape = 21, stroke = 0.8, color = "black") +
    scale_fill_gradientn(colours = c("lightgrey", "#4E68A3", "#DDB745", "#E67336", "#DA2F2F"), name = "RNA Exp (scaled)") +
    scale_size_continuous(range = c(2, 8), name = "-log10(P-value)") +
    labs(x = "Enriched Stage / Group", y = "Transcription Factor Motif") +
    theme_science() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# ------------------------------------------------------------------------------
# Example Execution
# ------------------------------------------------------------------------------
# target_tfs <- c("Nrf1", "Nr5a2", "Klf17", "Pou5f1", "Sox2", "Gata3", "Cdx2")
# p_integrated <- plot_ame_rna_integration("/path/to/ame_results", embryo.integrated3, "Stage", target_tfs)
# print(p_integrated)