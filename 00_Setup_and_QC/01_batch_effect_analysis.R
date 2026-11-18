# ==============================================================================
# Script: 01_batch_effect_analysis.R
# Purpose: Evaluate batch effects by comparing cell stage compositions across groups.
# ==============================================================================
source("00_utils.R")

# ------------------------------------------------------------------------------
# User Inputs (Adjust these for your specific project)
# ------------------------------------------------------------------------------
seurat_obj <- embryo.integrated3 # Replace with your merged Seurat object
group_col <- "Group"             # Column name for batch/group in meta.data
stage_col <- "Stage"             # Column name for developmental stage
stage_order <- opts$stage_lineage # Pre-defined stage order from 00_utils.R

# ------------------------------------------------------------------------------
# Data Preparation
# ------------------------------------------------------------------------------
meta <- seurat_obj@meta.data

# Ensure required columns exist
if(!(group_col %in% colnames(meta)) | !(stage_col %in% colnames(meta))) {
  stop("Specified group or stage columns not found in Seurat object metadata.")
}

# Factorize Stage for ordered plotting
meta[[stage_col]] <- factor(meta[[stage_col]], levels = stage_order)

# Calculate Cell Numbers and Proportions per Group/Stage
stage_prop <- meta %>%
  group_by(.data[[group_col]], .data[[stage_col]]) %>%
  summarise(Cell_Number = n(), .groups = "drop") %>%
  group_by(.data[[group_col]]) %>%
  mutate(Proportion = Cell_Number / sum(Cell_Number)) %>%
  ungroup()

print(head(stage_prop))

# ------------------------------------------------------------------------------
# Visualization
# ------------------------------------------------------------------------------

# Plot 1: Stacked Barplot (Proportions)
p_prop <- ggplot(stage_prop, aes(x = .data[[group_col]], y = Proportion, fill = .data[[stage_col]])) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Stage Composition by Batch",
       x = "Batch (Group)", y = "Proportion", fill = "Stage") +
  theme_science() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 2: Stacked Barplot (Absolute Counts)
p_count <- ggplot(stage_prop, aes(x = .data[[group_col]], y = Cell_Number, fill = .data[[stage_col]])) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(title = "Stage Cell Numbers by Batch",
       x = "Batch (Group)", y = "Cell Number", fill = "Stage") +
  theme_science() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ------------------------------------------------------------------------------
# Display & Save
# ------------------------------------------------------------------------------
print(p_prop)
print(p_count)

# ggsave(paste0(opts$outdir, "/batch_effect_proportion.pdf"), plot = p_prop, width = 7, height = 5)
# ggsave(paste0(opts$outdir, "/batch_effect_counts.pdf"), plot = p_count, width = 7, height = 5)