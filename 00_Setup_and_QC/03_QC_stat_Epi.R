# ==============================================================================
# Script: 02_QC_stat_Epi.R
# Purpose: Quality control and filtering for single-cell Methylation and Accessibility data.
# ==============================================================================
source("00_utils.R")

# ------------------------------------------------------------------------------
# Custom Plotting Theme Override
# ------------------------------------------------------------------------------
barplot_theme <- function() {
  theme_science() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
}

# ------------------------------------------------------------------------------
# Data Loading & Initialization
# ------------------------------------------------------------------------------
# Ensure you have set opts$basedir correctly in 00_utils.R before running
meta_file <- paste0(opts$basedir, "/sample_metadata.csv")
dir_met <- paste0(opts$basedir, "/met/cpg_level")
dir_acc <- paste0(opts$basedir, "/acc/cpg_level")

# Load structured metadata
md_info <- load_metadata(meta_file, id_columns = c("id_met", "id_acc"))
meta_df <- md_info$metadata
cells_met <- md_info$valid_cells$id_met
cells_acc <- md_info$valid_cells$id_acc

# Thresholds (imported from opts or define locally)
thr_met_cov <- 2.5e4
thr_acc_cov <- 5e5
thr_acc_mean_low <- 0.3
thr_acc_mean_high <- 0.8

# ------------------------------------------------------------------------------
# Metric Extraction Function
# ------------------------------------------------------------------------------
calculate_coverage <- function(cells, data_dir, ext) {
  stats <- data.table(id = cells, coverage = as.numeric(NA), mean_rate = as.numeric(NA))
  
  for (cell in cells) {
    file_path <- sprintf("%s/%s.%s.txt", data_dir, cell, ext)
    if (file.exists(file_path)) {
      tmp <- fread(file_path) %>% .[, c(1,2,5)] %>% setnames(c("chr", "pos", "rate"))
      stats[id == cell, `:=`(coverage = nrow(tmp), mean_rate = mean(tmp$rate, na.rm = TRUE))]
    }
  }
  return(stats[!is.na(coverage)])
}

# Calculate Stats
cat("Processing Methylation files...\n")
stat_met <- calculate_coverage(cells_met, dir_met, "met")

cat("Processing Accessibility files...\n")
stat_acc <- calculate_coverage(cells_acc, dir_acc, "acc")

# ------------------------------------------------------------------------------
# Quality Control Plotting & Filtering
# ------------------------------------------------------------------------------

# 1. Methylation QC
stat_met <- stat_met[order(coverage)]
stat_met[, id := factor(id, levels = id)]
fail_met <- stat_met[coverage < thr_met_cov, id]

p_met <- ggplot(stat_met, aes(x = id, y = coverage)) +
  geom_bar(stat = "identity", fill = "#F8766D", color = "#F8766D") +
  geom_hline(yintercept = thr_met_cov, colour = "black", linetype = "dashed") +
  labs(y = "Number of observed CpG sites") +
  scale_y_continuous(expand = c(0,0), limits = c(0, max(stat_met$coverage) * 1.1)) +
  barplot_theme()

# 2. Accessibility QC
stat_acc <- stat_acc[order(coverage)]
stat_acc[, id := factor(id, levels = id)]
fail_acc <- stat_acc[coverage < thr_acc_cov | mean_rate < thr_acc_mean_low | mean_rate > thr_acc_mean_high, id]

p_acc_cov <- ggplot(stat_acc, aes(x = id, y = coverage)) +
  geom_bar(stat = "identity", fill = "#00BFC4", color = "#00BFC4") +
  geom_hline(yintercept = thr_acc_cov, colour = "black", linetype = "dashed") +
  labs(y = "Number of observed GpC sites") +
  scale_y_continuous(expand = c(0,0), limits = c(0, max(stat_acc$coverage) * 1.1)) +
  barplot_theme()

p_acc_mean <- ggplot(stat_acc[order(mean_rate)], aes(x = id, y = mean_rate)) +
  geom_bar(stat = "identity", fill = "#00BFC4", color = "#00BFC4") +
  geom_hline(yintercept = c(thr_acc_mean_low, thr_acc_mean_high), colour = "black", linetype = "dashed") +
  labs(y = "Mean accessibility rate") +
  barplot_theme()

# ------------------------------------------------------------------------------
# Metadata Update & Output
# ------------------------------------------------------------------------------
print(p_met)
print(p_acc_cov)
print(p_acc_mean)

cat(sprintf("Cells failing Met QC: %d\nCells failing Acc QC: %d\n", length(fail_met), length(fail_acc)))

# Update global metadata table
meta_df[, pass_metQC := ifelse(id_met %in% fail_met, FALSE, TRUE)]
meta_df[, pass_accQC := ifelse(id_acc %in% fail_acc, FALSE, TRUE)]

# write.csv(stat_met, paste0(opts$outdir, "/sample_stats_met.csv"), row.names = FALSE)
# write.csv(stat_acc, paste0(opts$outdir, "/sample_stats_acc.csv"), row.names = FALSE)
# write.csv(meta_df, paste0(opts$outdir, "/sample_metadata_QC_filtered.csv"), row.names = FALSE)