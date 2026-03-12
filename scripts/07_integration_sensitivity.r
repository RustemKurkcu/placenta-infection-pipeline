source("R/config.R")
source("R/helpers_io.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(qs)
})

run_harmony_compat <- function(seu_obj, batch_col, dims = 1:30) {
  fn <- get("RunHarmony", envir = asNamespace("harmony"))
  fn_formals <- names(formals(fn))
  
  args <- list(
    object = seu_obj,
    group.by.vars = batch_col,
    verbose = FALSE
  )
  
  # Handle API differences between harmony versions.
  if ("reduction" %in% fn_formals) {
    args$reduction <- "pca"
  } else if ("reduction.use" %in% fn_formals) {
    args$reduction.use <- "pca"
  }
  
  if ("dims.use" %in% fn_formals) args$dims.use <- dims
  if ("assay.use" %in% fn_formals) args$assay.use <- "RNA"
  
  do.call(fn, args)
}

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("07_integration_sensitivity starting.", log_file = log_file)

seu <- qs::qread(file.path(CFG$dirs$objects, "seu_with_scores.qs"))
DefaultAssay(seu) <- "RNA"

# Use donor as batch variable by default; adjust in config if needed.
batch_col <- CFG$cols$donor_id
if (!batch_col %in% colnames(seu@meta.data)) stop("Missing donor/batch column: ", batch_col)

# Base preprocessing for PCA
seu <- NormalizeData(seu, verbose = FALSE)

# Guard against memory spikes on large objects.
# Base preprocessing for PCA
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures =4000, verbose = FALSE)
seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 50, verbose = FALSE)

# 1) Non-integrated reference embedding
seu <- RunUMAP(seu, dims = 1:30, reduction = "pca", reduction.name = "umap_rna", verbose = FALSE)

# 2) Harmony embedding (if package is available)
if (requireNamespace("harmony", quietly = TRUE)) {
  seu <- run_harmony_compat(seu_obj = seu, batch_col = batch_col, dims = 1:30)
  seu <- RunUMAP(seu, dims = 1:30, reduction = "harmony", reduction.name = "umap_harmony", verbose = FALSE)
  
  p_h1 <- DimPlot(seu, reduction = "umap_harmony", group.by = CFG$cols$cell_type, raster = TRUE) +
    ggtitle("Harmony UMAP by cell type")
  p_h2 <- DimPlot(seu, reduction = "umap_harmony", group.by = batch_col, raster = TRUE) +
    ggtitle("Harmony UMAP by donor")
  p_h3 <- DimPlot(seu, reduction = "umap_harmony", group.by = CFG$cols$infection, raster = TRUE) +
    ggtitle("Harmony UMAP by infection")
  
  save_plot(file.path(CFG$dirs$figures, "Fig16_Harmony_UMAP_cell_type.pdf"), p_h1, w = 12, h = 8)
  save_plot(file.path(CFG$dirs$figures, "Fig17_Harmony_UMAP_donor.pdf"), p_h2, w = 12, h = 8)
  save_plot(file.path(CFG$dirs$figures, "Fig18_Harmony_UMAP_infection.pdf"), p_h3, w = 12, h = 8)
  
  write_legend(
    "Fig16-18", "Integration sensitivity with Harmony",
    hypothesis = "Cross-donor comparability can improve after batch correction without erasing infection biology.",
    methods = "Run Harmony on PCA using donor as batch variable; compare cell-type, donor, and infection overlays.",
    readout = "Desired pattern: reduced donor-driven segregation while preserving cell-type structure and condition-associated shifts.",
    interpretation_template = "- If donor separation dominates after correction, revisit QC and batch variables.\n- If infection signal disappears globally, correction may be too aggressive.",
    outfile = file.path(CFG$dirs$legends, "Fig16_18_legend.md")
  )
} else {
  log_msg("Harmony not installed; skipping Harmony integration.", log_file = log_file)
}

# Save updated object with extra reductions for downstream inspection
qs::qsave(seu, file.path(CFG$dirs$objects, "seu_with_integration_checks.qs"))
log_msg("07_integration_sensitivity done.", log_file = log_file)