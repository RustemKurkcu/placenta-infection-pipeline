source("R/config.R")
source("R/helpers_io.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(qs)
  library(patchwork)
})

# =========================================================
# 07b_parameter_sensitivity.R
#
# PURPOSE
# -------
# Test a small, biologically reasonable grid of:
#   - HVGs: 3000, 4000, 5000
#   - PCs : 30, 40, 50
#
# Why this exists:
#   1) Placenta is highly heterogeneous, so 2000 HVGs may be too conservative.
#   2) We want an evidence-based choice of dimensionality for manuscript-quality work.
#   3) We want to preserve biologically meaningful structure:
#        - major placental cell types
#        - infection/time patterns
#        - donor correction without over-correction
#        - FIB2 / stromal / trophoblast / immune structure
#
# DESIGN CHOICES
# --------------
# - Keep ElbowPlot diagnostics
# - Skip JackStraw for now (computationally expensive)
# - Save intermediate objects with qs for reproducibility and recovery
# - Remove temporary objects and run gc() inside loops to reduce memory pressure
# - Use default Harmony reduction name "harmony" within each loop
#   because each loop starts from a fresh object
# =========================================================

# -----------------------------
# SETTINGS
# -----------------------------
log_file <- file.path(CFG$dirs$logs, "pipeline.log")
out_tbl  <- file.path(CFG$dirs$tables, "integration_param_sweep_summary.csv")
obj_dir  <- file.path(CFG$dirs$objects, "parameter_sweep")
dir.create(obj_dir, showWarnings = FALSE, recursive = TRUE)

# Reduced grid: practical + biologically sensible first pass
hvg_grid <- c(3000, 4000, 5000)
pc_grid  <- c(30, 40, 50)

seed <- 42

# Save per-combination objects?
# TRUE = very reproducible, but more disk usage
# FALSE = lighter on disk
save_combo_objects <- TRUE

# -----------------------------
# HELPERS
# -----------------------------

safe_dimplot <- function(seu_obj, reduction, group_by, title_txt) {
  if (!reduction %in% names(seu_obj@reductions)) return(NULL)
  if (!group_by %in% colnames(seu_obj@meta.data)) return(NULL)
  
  DimPlot(seu_obj, reduction = reduction, group.by = group_by, raster = TRUE) +
    ggtitle(title_txt)
}

save_if_not_null <- function(path, p, w = 10, h = 8) {
  if (!is.null(p)) save_plot(path, p, w = w, h = h)
}

# Simple elbow suggestion:
# returns the first PC after which marginal gain becomes small
# This is a heuristic, not a substitute for biological review.
suggest_elbow_pc <- function(seu_obj, min_pc = 10, max_pc = 50, rel_drop_thresh = 0.10) {
  stdev <- seu_obj[["pca"]]@stdev
  max_pc <- min(max_pc, length(stdev))
  
  if (max_pc < (min_pc + 2)) return(NA_integer_)
  
  vars <- stdev[1:max_pc]^2
  deltas <- -diff(vars)  # variance drop from PC_i to PC_{i+1}
  
  # relative drop compared with the first drop
  rel <- deltas / deltas[1]
  
  # first PC where the incremental drop becomes "small"
  idx <- which(seq_along(rel) >= min_pc & rel < rel_drop_thresh)[1]
  
  if (is.na(idx)) {
    return(max_pc)
  } else {
    return(idx)
  }
}

prep_pca <- function(obj, nfeatures, npcs = 50) {
  DefaultAssay(obj) <- "RNA"
  
  # Drop old reductions/graphs to keep the working copy lighter.
  # Seurat provides DietSeurat specifically for slimming objects. 
  obj <- DietSeurat(
    object = obj,
    assays = "RNA",
    layers = c("counts", "data"),
    dimreducs = NULL,
    graphs = NULL,
    misc = TRUE
  )
  
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = nfeatures,
    verbose = FALSE
  )
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = npcs, verbose = FALSE)
  
  obj
}

pc_var_table <- function(seu_obj) {
  stdev <- seu_obj[["pca"]]@stdev
  data.frame(
    PC = seq_along(stdev),
    stdev = stdev,
    var_explained = (stdev^2) / sum(stdev^2),
    cumvar = cumsum((stdev^2) / sum(stdev^2))
  )
}

run_harmony_safe <- function(seu_obj, batch_col, dims_use) {
  # Harmony documentation indicates PCA must already exist.
  # We use the default reduction.save = "harmony" within each fresh loop object.
  harmony::RunHarmony(
    object = seu_obj,
    group.by.vars = batch_col,
    reduction = "pca",
    dims.use = dims_use,
    assay.use = "RNA",
    reduction.save = "harmony",
    verbose = FALSE
  )
}

cleanup_big_objects <- function() {
  gc(verbose = FALSE)
}

# -----------------------------
# LOAD OBJECT
# -----------------------------
log_msg("07b_parameter_sensitivity starting.", log_file = log_file)

seu <- qs::qread(file.path(CFG$dirs$objects, "seu_with_scores.qs"))
DefaultAssay(seu) <- "RNA"

batch_col <- CFG$cols$donor_id
if (!batch_col %in% colnames(seu@meta.data)) {
  stop("Missing donor/batch column: ", batch_col)
}

# Build combined condition label if absent.
if (!"condition" %in% colnames(seu@meta.data) &&
    all(c(CFG$cols$infection, CFG$cols$hpi) %in% colnames(seu@meta.data))) {
  seu$condition <- paste0(seu[[CFG$cols$infection]][,1], "_", seu[[CFG$cols$hpi]][,1])
}

summary_rows <- list()

# Save the base object metadata snapshot for reproducibility
qs::qsave(
  x = seu,
  file = file.path(obj_dir, "base_input_object_for_parameter_sweep.qs")
)

# -----------------------------
# PARAMETER SWEEP
# -----------------------------
for (n_hvg in hvg_grid) {
  
  log_msg("Starting HVG setting: ", n_hvg, log_file = log_file)
  
  # Fresh working copy for this HVG setting
  obj <- tryCatch(
    prep_pca(seu, nfeatures = n_hvg, npcs = max(pc_grid)),
    error = function(e) {
      log_msg("prep_pca failed for HVG=", n_hvg, " :: ", conditionMessage(e), log_file = log_file)
      stop(e)
    }
  )
  
  pca_tbl <- pc_var_table(obj)
  elbow_pc_auto <- suggest_elbow_pc(obj, min_pc = 10, max_pc = max(pc_grid), rel_drop_thresh = 0.10)
  
  # Save PCA variance table
  pca_tbl$hvg <- n_hvg
  pca_tbl$auto_elbow_pc <- elbow_pc_auto
  
  write.csv(
    pca_tbl,
    file.path(CFG$dirs$tables, paste0("PCA_variance_HVG", n_hvg, ".csv")),
    row.names = FALSE
  )
  
  # Save HVG-level object after PCA
  qs::qsave(
    x = obj,
    file = file.path(obj_dir, paste0("seu_PCA_HVG", n_hvg, ".qs"))
  )
  
  # Elbow plot with automatic suggestion annotated
  p_elbow <- ElbowPlot(obj, ndims = max(pc_grid), reduction = "pca") +
    geom_vline(xintercept = elbow_pc_auto, linetype = "dashed") +
    ggtitle(paste0("ElbowPlot | HVGs=", n_hvg, " | auto elbow ~ PC ", elbow_pc_auto))
  
  save_plot(
    file.path(CFG$dirs$figures, paste0("FigParam_Elbow_HVG", n_hvg, ".pdf")),
    p_elbow, w = 9, h = 6
  )
  
  # Variable feature plot
  p_vf <- VariableFeaturePlot(obj) +
    ggtitle(paste0("VariableFeaturePlot | HVGs=", n_hvg))
  
  save_plot(
    file.path(CFG$dirs$figures, paste0("FigParam_VariableFeatures_HVG", n_hvg, ".pdf")),
    p_vf, w = 9, h = 6
  )
  
  # Per-PC loop
  for (n_pc in pc_grid) {
    
    log_msg("Running HVG=", n_hvg, " | PCs=", n_pc, log_file = log_file)
    
    # Fresh copy from the HVG-specific PCA object
    tmp <- obj
    
    # -----------------------------
    # Non-integrated UMAP
    # -----------------------------
    tmp <- RunUMAP(
      tmp,
      dims = 1:n_pc,
      reduction = "pca",
      reduction.name = "umap_rna",
      n.neighbors = 30,
      min.dist = 0.3,
      verbose = FALSE
    )
    
    # -----------------------------
    # Harmony + Harmony UMAP
    # -----------------------------
    if (requireNamespace("harmony", quietly = TRUE)) {
      tmp <- tryCatch(
        run_harmony_safe(tmp, batch_col = batch_col, dims_use = 1:n_pc),
        error = function(e) {
          log_msg("Harmony failed for HVG=", n_hvg, " | PC=", n_pc, " :: ", conditionMessage(e), log_file = log_file)
          tmp
        }
      )
      
      if ("harmony" %in% names(tmp@reductions)) {
        tmp <- RunUMAP(
          tmp,
          dims = 1:n_pc,
          reduction = "harmony",
          reduction.name = "umap_harmony",
          n.neighbors = 30,
          min.dist = 0.3,
          verbose = FALSE
        )
      }
    }
    
    # -----------------------------
    # Save object for this combination
    # -----------------------------
    if (isTRUE(save_combo_objects)) {
      qs::qsave(
        x = tmp,
        file = file.path(obj_dir, paste0("seu_HVG", n_hvg, "_PC", n_pc, ".qs"))
      )
    }
    
    # -----------------------------
    # RNA UMAP panel
    # -----------------------------
    p1 <- safe_dimplot(tmp, "umap_rna", CFG$cols$cell_type,
                       paste0("RNA UMAP | cell type | HVG=", n_hvg, " PC=", n_pc))
    p2 <- safe_dimplot(tmp, "umap_rna", batch_col,
                       paste0("RNA UMAP | donor | HVG=", n_hvg, " PC=", n_pc))
    p3 <- safe_dimplot(tmp, "umap_rna", CFG$cols$infection,
                       paste0("RNA UMAP | infection | HVG=", n_hvg, " PC=", n_pc))
    p4 <- safe_dimplot(tmp, "umap_rna", "condition",
                       paste0("RNA UMAP | condition | HVG=", n_hvg, " PC=", n_pc))
    
    if (!is.null(p1) && !is.null(p2) && !is.null(p3) && !is.null(p4)) {
      panel <- (p1 | p2) / (p3 | p4)
      save_plot(
        file.path(CFG$dirs$figures, paste0("FigParam_RNA_HVG", n_hvg, "_PC", n_pc, ".pdf")),
        panel, w = 16, h = 12
      )
    }
    
    # -----------------------------
    # Harmony UMAP panel
    # -----------------------------
    if ("umap_harmony" %in% names(tmp@reductions)) {
      h1 <- safe_dimplot(tmp, "umap_harmony", CFG$cols$cell_type,
                         paste0("Harmony UMAP | cell type | HVG=", n_hvg, " PC=", n_pc))
      h2 <- safe_dimplot(tmp, "umap_harmony", batch_col,
                         paste0("Harmony UMAP | donor | HVG=", n_hvg, " PC=", n_pc))
      h3 <- safe_dimplot(tmp, "umap_harmony", CFG$cols$infection,
                         paste0("Harmony UMAP | infection | HVG=", n_hvg, " PC=", n_pc))
      h4 <- safe_dimplot(tmp, "umap_harmony", "condition",
                         paste0("Harmony UMAP | condition | HVG=", n_hvg, " PC=", n_pc))
      
      if (!is.null(h1) && !is.null(h2) && !is.null(h3) && !is.null(h4)) {
        hpanel <- (h1 | h2) / (h3 | h4)
        save_plot(
          file.path(CFG$dirs$figures, paste0("FigParam_Harmony_HVG", n_hvg, "_PC", n_pc, ".pdf")),
          hpanel, w = 16, h = 12
        )
      }
    }
    
    # -----------------------------
    # Summary row
    # -----------------------------
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      n_cells = ncol(tmp),
      hvg = n_hvg,
      pcs = n_pc,
      auto_elbow_pc = elbow_pc_auto,
      pca_var_pc1_10 = sum(pca_tbl$var_explained[1:min(10, nrow(pca_tbl))]),
      pca_var_pc1_20 = sum(pca_tbl$var_explained[1:min(20, nrow(pca_tbl))]),
      pca_var_pc1_n  = sum(pca_tbl$var_explained[1:min(n_pc, nrow(pca_tbl))]),
      harmony_success = "harmony" %in% names(tmp@reductions),
      stringsAsFactors = FALSE
    )
    
    # -----------------------------
    # Clean up large temp objects
    # -----------------------------
    rm(list = intersect(
      c("tmp", "p1", "p2", "p3", "p4", "panel", "h1", "h2", "h3", "h4", "hpanel"),
      ls()
    ))
    cleanup_big_objects()
  }
  
  # Clean HVG-level objects before next HVG
  rm(list = intersect(c("obj", "pca_tbl", "p_elbow", "p_vf", "elbow_pc_auto"), ls()))
  cleanup_big_objects()
}

# -----------------------------
# WRITE SUMMARY TABLE
# -----------------------------
summary_df <- bind_rows(summary_rows)
write.csv(summary_df, out_tbl, row.names = FALSE)

# Also save summary as qs
qs::qsave(summary_df, file.path(obj_dir, "integration_param_sweep_summary.qs"))

# -----------------------------
# SUMMARY PLOTS
# -----------------------------
p_var <- ggplot(summary_df, aes(x = pcs, y = pca_var_pc1_n, color = factor(hvg), group = hvg)) +
  geom_line() + geom_point() +
  labs(
    title = "Cumulative PCA variance explained across HVG/PC settings",
    x = "PCs used",
    y = "Cumulative variance explained",
    color = "HVGs"
  ) +
  theme_bw()

p_elbow_suggest <- ggplot(summary_df, aes(x = factor(hvg), y = auto_elbow_pc)) +
  geom_point(size = 3) +
  labs(
    title = "Automatic elbow suggestion by HVG setting",
    x = "HVGs",
    y = "Suggested elbow PC"
  ) +
  theme_bw()

save_plot(file.path(CFG$dirs$figures, "FigParam_Summary_variance_explained.pdf"), p_var, w = 10, h = 7)
save_plot(file.path(CFG$dirs$figures, "FigParam_Summary_auto_elbow.pdf"), p_elbow_suggest, w = 8, h = 6)

# -----------------------------
# LEGEND
# -----------------------------
write_legend(
  "FigParam", "Parameter sensitivity for HVGs and PCs",
  hypothesis = paste(
    "The selected dimensionality should preserve major placental lineages, infection/time-associated biology,",
    "and donor-corrected structure without evidence of over-correction."
  ),
  methods = paste(
    "Tested HVGs =", paste(hvg_grid, collapse = ","),
    "; tested PCs =", paste(pc_grid, collapse = ","),
    ". For each HVG setting, PCA was computed, ElbowPlot and VariableFeaturePlot were saved,",
    "and non-integrated plus Harmony-corrected UMAPs were generated for each PC setting.",
    "Intermediate objects were saved with qs for reproducibility and recovery."
  ),
  readout = paste(
    "Preferred settings should preserve biologically expected structure across cell type, donor, infection, and condition overlays.",
    "Automatic elbow suggestions are heuristic only and should be interpreted together with the UMAP panels and biological coherence."
  ),
  interpretation_template = paste(
    "- Prefer the smallest PC setting that preserves stable biological structure.",
    "\n- Favor HVG settings that preserve rare or important populations without obvious fragmentation.",
    "\n- If donor correction improves but infection/time structure collapses, Harmony may be over-correcting.",
    "\n- If 4000 and 5000 perform similarly, the simpler or more stable setting may be preferable."
  ),
  outfile = file.path(CFG$dirs$legends, "FigParam_legend.md")
)

# -----------------------------
# OPTIONAL JACKSTRAW (DEFERRED)
# -----------------------------
# Intentionally skipped for now.
# Seurat documentation notes that ElbowPlot is much faster than JackStraw
# and often corresponds well with significant dimensions.
# If needed, run JackStraw later only on 1-2 finalist settings and on a subset.

log_msg("07b_parameter_sensitivity done.", log_file = log_file)

