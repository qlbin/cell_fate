# ==============================================================================
# Script: SingleCell_Epi_Differential.R
# Purpose: Single-cell differential accessibility/methylation (Binomial/T-test) 
#          and Bulk DiffBind analysis integration.
# ==============================================================================
source("00_utils.R")
suppressPackageStartupMessages({
  library(data.table)
  library(purrr)
  library(DiffBind)
})

# ------------------------------------------------------------------------------
# 1. Single-Cell Feature Differential Analysis
# ------------------------------------------------------------------------------
run_sc_epi_differential <- function(data_dir, metadata_df, id_col, annos_list, 
                                    groupA, groupB, group_col = "Celltype10",
                                    test_method = "binomial", 
                                    min_cells = 2, min_cov = 5, fdr_threshold = 0.05, diff_threshold = 1) {
  
  # Extract relevant cells
  valid_cells <- metadata_df %>% filter(!!sym(group_col) %in% c(groupA, groupB)) %>% pull(!!sym(id_col))
  
  # Load and bind requested annotations
  data <- lapply(annos_list, function(n) {
    fread(cmd = sprintf("zcat < %s/%s.tsv.gz", data_dir, n), showProgress = FALSE, header = FALSE)
  }) %>% rbindlist() %>% setnames(c("id_cell", "id", "anno", "N_hit", "N_total", "rate"))
  
  # Merge and assign groups
  data <- data[id_cell %in% valid_cells]
  meta_subset <- metadata_df[metadata_df[[id_col]] %in% valid_cells, c(id_col, group_col), with = FALSE]
  setnames(meta_subset, id_col, "id_cell")
  
  data <- merge(data, meta_subset, by = "id_cell")
  data[, group := ifelse(get(group_col) == groupA, "A", "B")]
  
  # Coverage Filtering
  data <- data[N_total >= min_cov]
  
  # Ensure feature is detected in both groups
  data <- data[, Ngroup := length(unique(group)), by = .(id, anno)][Ngroup == 2][, Ngroup := NULL]
  
  # Minimum Cell Filtering
  remove_sites <- data[, .(N = min(.N)), by = .(id, group)][N < min_cells, unique(id)]
  data <- data[!id %in% remove_sites]
  
  # Statistical Testing
  if (test_method == "binomial") {
    diff_res <- data[, .(
      A_hit = sum(.SD[group == "A", N_hit]), A_miss = sum(.SD[group == "A", N_total - N_hit]),
      B_hit = sum(.SD[group == "B", N_hit]), B_miss = sum(.SD[group == "B", N_total - N_hit])
    ), by = .(id, anno)]
    
    diff_res[, p.value := fisher.test(matrix(c(A_hit, A_miss, B_hit, B_miss), 2, 2))[["p.value"]], by = .(id, anno)]
    diff_res[, `:=`(rateA = 100 * (A_hit / (A_hit + A_miss)), rateB = 100 * (B_hit / (B_hit + B_miss)))]
  } 
  
  # Multiple Testing Correction
  diff_res[, diff := rateB - rateA]
  diff_res[, padj_fdr := p.adjust(p.value, method = "fdr"), by = anno]
  diff_res[, sig := (p.value <= fdr_threshold & abs(diff) > diff_threshold)]
  diff_res <- diff_res[order(padj_fdr)]
  
  return(diff_res)
}

# ------------------------------------------------------------------------------
# 2. Bulk DiffBind Routine Wrapper
# ------------------------------------------------------------------------------
run_bulk_diffbind <- function(sample_sheet_csv, contrast_col = DBA_CONDITION, min_members = 2) {
  
  samples <- read.csv(sample_sheet_csv)
  dbObj <- dba(sampleSheet = samples)
  dbObj <- dba.count(dbObj, bUseSummarizeOverlaps = TRUE)
  dbObj <- dba.normalize(dbObj)
  
  dbObj <- dba.contrast(dbObj, minMembers = min_members, categories = contrast_col)
  dbObj <- dba.analyze(dbObj, method = DBA_DESEQ2)
  
  return(dbObj)
}