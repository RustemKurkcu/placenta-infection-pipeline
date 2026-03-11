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

downsample_for_speed <- function(seu_obj, max_cells = 80000, seed = 42) {
  n <- ncol(seu_obj)
  if (n <= max_cells) return(seu_obj)

  set.seed(seed)
  md <- seu_obj@meta.data
  key <- if (all(c(CFG$cols$cell_type, CFG$cols$donor_id) %in% colnames(md))) {
    paste(md[[CFG$cols$cell_type]], md[[CFG$cols$donor_id]], sep = "__")
  } else {
    rep("all", n)
  }

  split_idx <- split(seq_len(n), key)
  per_group <- max(1L, floor(max_cells / length(split_idx)))
  keep <- unlist(lapply(split_idx, function(ix) sample(ix, size = min(length(ix), per_group))), use.names = FALSE)

  if (length(keep) > max_cells) keep <- sample(keep, size = max_cells)
  subset(seu_obj, cells = colnames(seu_obj)[sort(keep)])
}

safe_dimplot <- function(seu_obj, reduction, group_by, title_txt) {
  if (!group_by %in% colnames(seu_obj@meta.data)) return(NULL)
  DimPlot(seu_obj, reduction = reduction, group.by = group_by, raster = TRUE) +
    ggtitle(title_txt)
}

save_if_not_null <- function(path, p, w = 12, h = 8) {
  if (!is.null(p)) save_plot(path, p, w = w, h = h)
}

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("07_integration_sensitivity starting.", log_file = log_file)

seu <- qs::qread(file.path(CFG$dirs$objects, "seu_with_scores.qs"))
DefaultAssay(seu) <- "RNA"

batch_col <- CFG$cols$donor_id
if (!batch_col %in% colnames(seu@meta.data)) stop("Missing donor/batch column: ", batch_col)

# Create condition if missing.
if (!"condition" %in% colnames(seu@meta.data) &&
    all(c(CFG$cols$infection, CFG$cols$hpi) %in% colnames(seu@meta.data))) {
  seu$condition <- paste0(seu[[CFG$cols$infection]][,1], "_", seu[[CFG$cols$hpi]][,1])
}

seu_int <- downsample_for_speed(seu, max_cells = 80000, seed = 42)
log_msg("Integration sensitivity running on ", ncol(seu_int), " cells (from ", ncol(seu), ").", log_file = log_file)

# Base preprocessing for PCA
seu <- NormalizeData(seu, verbose = FALSE)

# Guard against memory spikes on large objects.
# Base preprocessing for PCA
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures =4000, verbose = FALSE)
seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 50, verbose = FALSE)

# 1) Non-integrated embeddings
seu_int <- RunUMAP(
  seu_int,
  dims = 1:30,
  reduction = "pca",
  reduction.name = "umap_rna",
  n.neighbors = 30,
  min.dist = 0.3,
  verbose = FALSE
)

if (requireNamespace("Rtsne", quietly = TRUE)) {
  seu_int <- RunTSNE(seu_int, dims = 1:30, reduction = "pca", reduction.name = "tsne_rna", check_duplicates = FALSE, verbose = FALSE)
} else {
  log_msg("Rtsne not installed; skipping tSNE on RNA PCA.", log_file = log_file)
}

# 2) Harmony embeddings
if (requireNamespace("harmony", quietly = TRUE)) {
  seu_int <- withCallingHandlers(
    run_harmony_compat(seu_obj = seu_int, batch_col = batch_col, dims = 1:30),
    warning = function(w) {
      log_msg("Harmony warning: ", conditionMessage(w), log_file = log_file)
      invokeRestart("muffleWarning")
    }
  )

  seu_int <- RunUMAP(
    seu_int,
    dims = 1:30,
    reduction = "harmony",
    reduction.name = "umap_harmony",
    n.neighbors = 30,
    min.dist = 0.3,
    verbose = FALSE
  )

  if (requireNamespace("Rtsne", quietly = TRUE)) {
    seu_int <- RunTSNE(seu_int, dims = 1:30, reduction = "harmony", reduction.name = "tsne_harmony", check_duplicates = FALSE, verbose = FALSE)
  }

  # Core Harmony plots
  save_if_not_null(file.path(CFG$dirs$figures, "Fig16_Harmony_UMAP_cell_type.pdf"),
                   safe_dimplot(seu_int, "umap_harmony", CFG$cols$cell_type, "Harmony UMAP by cell type"))
  save_if_not_null(file.path(CFG$dirs$figures, "Fig17_Harmony_UMAP_donor.pdf"),
                   safe_dimplot(seu_int, "umap_harmony", batch_col, "Harmony UMAP by donor"))
  save_if_not_null(file.path(CFG$dirs$figures, "Fig18_Harmony_UMAP_infection.pdf"),
                   safe_dimplot(seu_int, "umap_harmony", CFG$cols$infection, "Harmony UMAP by infection"))
  save_if_not_null(file.path(CFG$dirs$figures, "Fig19_Harmony_UMAP_stage.pdf"),
                   safe_dimplot(seu_int, "umap_harmony", CFG$cols$stage, "Harmony UMAP by stage"))
  save_if_not_null(file.path(CFG$dirs$figures, "Fig20_Harmony_UMAP_condition.pdf"),
                   safe_dimplot(seu_int, "umap_harmony", "condition", "Harmony UMAP by condition"))

  # tSNE companion plots
  if ("tsne_harmony" %in% names(seu_int@reductions)) {
    save_if_not_null(file.path(CFG$dirs$figures, "Fig21_Harmony_tSNE_cell_type.pdf"),
                     safe_dimplot(seu_int, "tsne_harmony", CFG$cols$cell_type, "Harmony tSNE by cell type"))
    save_if_not_null(file.path(CFG$dirs$figures, "Fig22_Harmony_tSNE_infection.pdf"),
                     safe_dimplot(seu_int, "tsne_harmony", CFG$cols$infection, "Harmony tSNE by infection"))
    save_if_not_null(file.path(CFG$dirs$figures, "Fig23_Harmony_tSNE_stage.pdf"),
                     safe_dimplot(seu_int, "tsne_harmony", CFG$cols$stage, "Harmony tSNE by stage"))
  }

  write_legend(
    "Fig16-23", "Integration sensitivity with Harmony + UMAP/tSNE",
    hypothesis = "Cross-donor comparability can improve after correction without erasing infection/time biology.",
    methods = "Downsample (if needed) for speed, run PCA, Harmony (donor batch), then UMAP and tSNE on Harmony embeddings; color by cell type, donor, infection, stage, and condition.",
    readout = "Desired: donor mixing improves while biological structure (cell type + infection/time patterns) remains interpretable. Compare UMAP (global manifold) with tSNE (local neighborhoods).",
    interpretation_template = "- UMAP is typically better for global structure trends; tSNE emphasizes local neighborhoods.\n- If donor still dominates, revisit QC/batch definition.\n- If infection/stage structure vanishes, integration may be over-correcting.",
    outfile = file.path(CFG$dirs$legends, "Fig16_23_legend.md")
  )
} else {
  log_msg("Harmony not installed; skipping Harmony integration.", log_file = log_file)
}

# Save updated object with extra reductions for downstream inspection
qs::qsave(seu_int, file.path(CFG$dirs$objects, "seu_with_integration_checks.qs"))
log_msg("07_integration_sensitivity done.", log_file = log_file)
