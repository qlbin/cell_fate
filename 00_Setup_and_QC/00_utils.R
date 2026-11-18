# ==============================================================================
# Script: 00_utils.R
# Purpose: Global configurations, standard libraries, and helper functions 
#          for multi-omics analysis (General template).
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Standard Libraries
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  # Core Data Manipulation
  library(data.table)
  library(purrr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  
  # Visualization
  library(ggplot2)
  library(ggpubr)
  library(ComplexHeatmap)
})

# ------------------------------------------------------------------------------
# 2. Publication-Ready ggplot2 Theme (Science Advances standard)
# ------------------------------------------------------------------------------
theme_science <- function(base_size = 12, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(size = rel(1.2), face = "bold", hjust = 0.5),
      axis.title = element_text(size = rel(1.1), face = "bold"),
      axis.text = element_text(size = rel(1.0), color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = rel(0.9)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.8),
      strip.text = element_text(face = "bold", size = rel(1.1))
    )
}

# ------------------------------------------------------------------------------
# 3. Global Options & Filtering Thresholds
# ------------------------------------------------------------------------------
opts <- list()

# Define Base I/O Directories (Modify before running project)
opts$basedir <- "./results"
opts$outdir <- paste0(opts$basedir, "/output")

# General filtering thresholds (Default defaults for sc-omics)
opts$met_min.CpGs <- 1        
opts$met_min.cells <- 15      
opts$met_nfeatures <- 15000   

opts$acc_min.GpCs <- 5        
opts$acc_min.cells <- 15      
opts$acc_nfeatures <- 15000   

opts$rna_min.cdr <- 0.25      
opts$rna_ngenes <- 4000       

# ------------------------------------------------------------------------------
# 4. Helper Functions
# ------------------------------------------------------------------------------

#' Load and format basic metadata
#' @param metadata_path Path to the metadata CSV
#' @param id_columns Character vector of column names representing cell IDs for different modalities
#' @return A list containing the full metadata and valid cell lists per modality
load_metadata <- function(metadata_path, id_columns = c("id_met", "id_acc", "id_rna")) {
  
  if(!file.exists(metadata_path)) stop("Metadata file not found: ", metadata_path)
  md <- fread(metadata_path)
  
  # Basic cleaning
  md <- na.omit(md)
  
  # Extract valid cells for each modality dynamically
  cells <- list()
  for(col in id_columns) {
    if(col %in% colnames(md)) {
      # Keep cells where the ID column is not NA
      cells[[col]] <- md[!is.na(get(col)), get(col)]
    } else {
      warning(paste("Column", col, "not found in metadata."))
    }
  }
  
  return(list(metadata = md, valid_cells = cells))
}

message("00_utils.R loaded successfully. Global options, generic functions, and theme_science() are available.")