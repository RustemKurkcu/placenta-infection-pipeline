source("R/config.R")
source("R/helpers_io.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(qs)
})

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("14_make_readable_embedding_figures starting.", log_file = log_file)

obj_candidates <- c(
  file.path(CFG$dirs$objects, "seu_with_architecture_transfer.qs"),
  file.path(CFG$dirs$objects, "seu_with_integration_checks.qs"),
  file.path(CFG$dirs$objects, "seu_with_scores.qs"),
  file.path(CFG$dirs$objects, "seu_clean.qs")
)
obj_path <- obj_candidates[file.exists(obj_candidates)][1]
if (is.na(obj_path) || length(obj_path) == 0) stop("No Seurat object found for script 14")

log_msg("Using object: ", obj_path, log_file = log_file)
seu <- qread(obj_path)

pick_reduction <- function(s) {
  reds <- names(s@reductions)
  if ("umap_harmony" %in% reds) return("umap_harmony")
  if ("umap_rna" %in% reds) return("umap_rna")
  if ("X_umap" %in% reds) return("X_umap")
  if ("umap" %in% reds) return("umap")
  NULL
}

red <- pick_reduction(seu)
if (is.null(red)) stop("No UMAP-like reduction found in object")

if (!"condition" %in% colnames(seu@meta.data) && all(c(CFG$cols$infection, CFG$cols$hpi) %in% colnames(seu@meta.data))) {
  seu$condition <- paste0(seu[[CFG$cols$infection]][,1], "_", seu[[CFG$cols$hpi]][,1])
}

plot_cols <- c(CFG$cols$cell_type, CFG$cols$infection, CFG$cols$hpi, CFG$cols$stage, CFG$cols$donor_id, "condition")
plot_cols <- unique(plot_cols[plot_cols %in% colnames(seu@meta.data)])

base_theme <- theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

for (grp in plot_cols) {
  p <- DimPlot(
    seu,
    reduction = red,
    group.by = grp,
    raster = TRUE,
    label = grp == CFG$cols$cell_type,
    repel = TRUE,
    pt.size = CFG$plotting$pt_size
  ) +
    ggtitle(paste0("Readable embedding by ", grp, " (", red, ")")) +
    base_theme

  save_plot(file.path(CFG$dirs$figures, paste0("FigReadable_", grp, "_", red, ".pdf")), p, w = 13, h = 9)
  save_plot(file.path(CFG$dirs$figures, paste0("FigReadable_", grp, "_", red, ".png")), p, w = 13, h = 9, dpi = 300)
}

if (all(c(CFG$cols$infection, CFG$cols$hpi) %in% colnames(seu@meta.data))) {
  grp <- CFG$cols$infection
  split_by <- CFG$cols$hpi
  p_split <- DimPlot(
    seu,
    reduction = red,
    group.by = grp,
    split.by = split_by,
    raster = TRUE,
    pt.size = CFG$plotting$pt_size,
    ncol = 3
  ) +
    ggtitle(paste0("Readable embedding: ", grp, " split by ", split_by, " (", red, ")")) +
    base_theme

  save_plot(file.path(CFG$dirs$figures, paste0("FigReadable_", grp, "_split_", split_by, "_", red, ".pdf")), p_split, w = 16, h = 9)
  save_plot(file.path(CFG$dirs$figures, paste0("FigReadable_", grp, "_split_", split_by, "_", red, ".png")), p_split, w = 16, h = 9, dpi = 300)
}

write_legend(
  "FigReadable", "Readable embedding panel set",
  hypothesis = "Improved formatting (larger canvas, clearer theme, consistent legends) increases interpretability of categorical embedding overlays.",
  methods = "Load best available Seurat object, choose best reduction, and regenerate key category plots with publication-readable sizing and theme settings.",
  readout = "Plots should have legible legends/text and enough panel size to distinguish category structure.",
  interpretation_template = "- Use these panels for human interpretation and manuscript assembly.\n- If overplotting remains high, facet or downsample by category.",
  outfile = file.path(CFG$dirs$legends, "FigReadable_legend.md")
)

log_msg("14_make_readable_embedding_figures done.", log_file = log_file)
