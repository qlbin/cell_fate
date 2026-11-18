# ==============================================================================
# Script: 01_convert_h5ad_to_Seurat.R
# Purpose: Convert Python AnnData exports (HDF5/TSV) into a processed Seurat object.
# ==============================================================================
source("00_utils.R")
suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
  library(Seurat)
})

# ------------------------------------------------------------------------------
# Data Construction Function
# ------------------------------------------------------------------------------
build_seurat_from_hdf5 <- function(h5_file, meta_file, umap_file, assay_name = "RNA") {
  
  if(!all(file.exists(c(h5_file, meta_file, umap_file)))) stop("Input files missing.")
  
  # Load HDF5 matrix data
  mydata <- h5read(h5_file, "mat")
  mat <- mydata$block0_values
  rownames(mat) <- mydata$axis0
  colnames(mat) <- mydata$axis1
  mat <- Matrix(mat, sparse = TRUE)
  
  # Load metadata and coords
  meta <- read.table(meta_file, sep="\t", header=TRUE, row.names=1)
  pos <- read.table(umap_file, sep="\t", header=TRUE, row.names=1)
  
  # Initialize
  obj <- CreateSeuratObject(counts = mat, assay = assay_name, meta.data = meta)
  
  # Temporary list to pass coords to clustering function
  obj@misc$precomputed_umap <- pos
  return(obj)
}

# ------------------------------------------------------------------------------
# Standard QC & Clustering Function
# ------------------------------------------------------------------------------
run_standard_seurat <- function(obj, mt_pattern = "^MT-", mt_list = NULL, 
                                dim_use = 30, res = 1,
                                nf_low = 200, nf_high = 6000, mt_cutoff = 10) {
  
  # MT annotation
  if (is.null(mt_list)) {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  } else {
    mt_list_valid <- intersect(mt_list, rownames(obj))
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, features = mt_list_valid)
  }
  
  # Filtering
  obj <- subset(obj, subset = nFeature_RNA > nf_low & 
                  percent.mt < mt_cutoff & nFeature_RNA < nf_high)
  
  # Processing Pipeline
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
  obj <- ScaleData(obj, features = rownames(obj))
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 40, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:dim_use)
  obj <- FindClusters(obj, resolution = res)
  obj <- RunUMAP(obj, dims = 1:dim_use)
  
  # Inject precomputed UMAP coords if available (optional replacement)
  if(!is.null(obj@misc$precomputed_umap)) {
    pos <- obj@misc$precomputed_umap
    common_cells <- intersect(rownames(pos), colnames(obj))
    if(length(common_cells) > 0) {
      obj@reductions$umap@cell.embeddings[common_cells, 1] <- pos[common_cells, 1]
      obj@reductions$umap@cell.embeddings[common_cells, 2] <- pos[common_cells, 2]
    }
  }
  
  return(obj)
}

# ------------------------------------------------------------------------------
# Example Usage 
# ------------------------------------------------------------------------------
# obj <- build_seurat_from_hdf5("mat.h5", "metadata.tsv", "position_X_umap.tsv", "RNA")
# obj <- run_standard_seurat(obj)
# saveRDS(obj, paste0(opts$outdir, "/Processed_Seurat.rds"))