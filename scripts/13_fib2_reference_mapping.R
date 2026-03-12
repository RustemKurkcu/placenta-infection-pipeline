source("R/config.R")
source("R/helpers_io.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(qs)
  library(Matrix)
})

# ============================================================
# 13_fib2_reference_mapping.R
#
# HYPOTHESIS
# ----------
# The broad stromal classes in the infection dataset (F, F_p, F_sm, PV)
# may contain a FIB2-like subpopulation similar to the architecture atlas.
#
# WHY THIS IMPROVES THE THESIS
# ----------------------------
# Your central 'metabolic reservoir' claim is strongest if you can map a
# FIB2-like compartment rather than only broad fibroblast labels.
#
# EXPECTED RESULT
# ---------------
# A subset of infection-dataset stromal / perivascular cells should map to
# FIB2-like labels with nontrivial prediction scores. If not, that weakens a
# strict FIB2 claim and suggests broader stromal states or missing reference
# resolution.
# ============================================================

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("13_fib2_reference_mapping starting.", log_file = log_file)

fig_dir <- CFG$dirs$figures
obj_dir <- CFG$dirs$objects
tbl_dir <- CFG$dirs$tables

# -----------------------------
# Helper functions
# -----------------------------
pick_plot_reduction <- function(seu_obj) {
  reds <- names(seu_obj@reductions)
  if ("umap_harmony" %in% reds) return("umap_harmony")
  if ("umap_rna" %in% reds) return("umap_rna")
  if ("X_umap" %in% reds) return("X_umap")
  if ("umap" %in% reds) return("umap")
  return(NULL)
}

prep_for_transfer <- function(obj, npcs = 30, nfeatures = 3000) {
  DefaultAssay(obj) <- "RNA"
  if (!"pca" %in% names(obj@reductions)) {
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = nfeatures, verbose = FALSE)
    obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
    obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = npcs, verbose = FALSE)
  }
  obj
}

guess_label_col <- function(md) {
  cand <- names(md)[vapply(md, function(x) is.character(x) || is.factor(x), logical(1))]
  fib_hits <- cand[
    vapply(cand, function(cc) {
      vals <- unique(as.character(md[[cc]]))
      any(grepl("FIB2|FIB1|Mat\\.FIB|fib", vals, ignore.case = TRUE))
    }, logical(1))
  ]
  if (length(fib_hits) > 0) return(fib_hits[1])
  stop("Could not identify a reference label column containing FIB-like labels.")
}

# -----------------------------
# Load infection object
# -----------------------------
seu_inf <- qs::qread(file.path(obj_dir, "seu_with_scores.qs"))
DefaultAssay(seu_inf) <- "RNA"

ct_col <- CFG$cols$cell_type

# Programmatic Seurat subsetting should use cells= rather than tidy-eval in subset=
keep_inf <- c("F", "F_p", "F_sm", "PV", "Endo_f", "EVT_1", "EVT_2", "iEVT", "VCT", "VCT_fusing")
keep_inf <- intersect(keep_inf, unique(as.character(seu_inf@meta.data[[ct_col]])))

cells_keep_inf <- rownames(seu_inf@meta.data)[
  as.character(seu_inf@meta.data[[ct_col]]) %in% keep_inf
]

if (length(cells_keep_inf) == 0) {
  stop("No infection cells matched keep_inf labels.")
}

seu_inf_sub <- subset(seu_inf, cells = cells_keep_inf)
seu_inf_sub <- prep_for_transfer(seu_inf_sub, npcs = 30, nfeatures = 3000)

# -----------------------------
# Locate architecture reference object
# EDIT THIS if your reference lives elsewhere
# -----------------------------
ref_candidates <- c(
  file.path("data", "architecture_reference.qs"),
  file.path("outputs", "objects", "architecture_reference.qs"),
  file.path("..", "hPlacenta-architecture", "outputs", "objects", "architecture_reference.qs"),
  file.path("..", "hPlacenta-architecture", "outputs", "objects", "seu_architecture.qs"),
  file.path("..", "hPlacenta-architecture", "outputs", "objects", "seu_hplacenta_architecture.qs")
)

ref_path <- ref_candidates[file.exists(ref_candidates)][1]

if (is.na(ref_path) || length(ref_path) == 0) {
  stop(
    "Could not find an architecture reference object.\n",
    "Please set ref_path manually to a .qs reference from the hPlacenta-architecture project."
  )
}

log_msg("Using architecture reference: ", ref_path, log_file = log_file)

seu_ref <- qs::qread(ref_path)
DefaultAssay(seu_ref) <- "RNA"

ref_label_col <- guess_label_col(seu_ref@meta.data)
log_msg("Using reference label column: ", ref_label_col, log_file = log_file)

# Keep a broad neighborhood around FIB2-like states
ref_vals <- unique(as.character(seu_ref@meta.data[[ref_label_col]]))
keep_ref <- ref_vals[
  grepl("FIB|fib|strom|peri|endo|EVT|VCT", ref_vals, ignore.case = TRUE)
]

if (length(keep_ref) == 0) {
  keep_ref <- ref_vals
}

cells_keep_ref <- rownames(seu_ref@meta.data)[
  as.character(seu_ref@meta.data[[ref_label_col]]) %in% keep_ref
]

seu_ref_sub <- subset(seu_ref, cells = cells_keep_ref)
seu_ref_sub <- prep_for_transfer(seu_ref_sub, npcs = 30, nfeatures = 3000)

# -----------------------------
# Transfer labels
# -----------------------------
anchors <- FindTransferAnchors(
  reference = seu_ref_sub,
  query = seu_inf_sub,
  normalization.method = "LogNormalize",
  reduction = "pcaproject",
  dims = 1:30,
  verbose = FALSE
)

pred <- TransferData(
  anchorset = anchors,
  refdata = seu_ref_sub@meta.data[[ref_label_col]],
  dims = 1:30,
  verbose = FALSE
)

# pred typically contains predicted.id and prediction.score.max
colnames(pred) <- sub("^predicted.id$", "arch_predicted_label", colnames(pred))
colnames(pred) <- sub("^prediction.score.max$", "arch_prediction_score_max", colnames(pred))

seu_inf_sub <- AddMetaData(seu_inf_sub, metadata = pred)

# Map predictions back onto full infection object
seu_inf$arch_predicted_label <- NA_character_
seu_inf$arch_prediction_score_max <- NA_real_

common_cells <- intersect(colnames(seu_inf), colnames(seu_inf_sub))
seu_inf$arch_predicted_label[match(common_cells, colnames(seu_inf))] <-
  seu_inf_sub$arch_predicted_label[match(common_cells, colnames(seu_inf_sub))]
seu_inf$arch_prediction_score_max[match(common_cells, colnames(seu_inf))] <-
  seu_inf_sub$arch_prediction_score_max[match(common_cells, colnames(seu_inf_sub))]

# Save mapped object
qs::qsave(seu_inf, file.path(obj_dir, "seu_with_architecture_transfer.qs"), preset = "high")

# -----------------------------
# Summary table
# -----------------------------
summary_tbl <- seu_inf@meta.data %>%
  filter(!is.na(arch_predicted_label)) %>%
  count(.data[[ct_col]], arch_predicted_label, name = "n_cells") %>%
  arrange(desc(n_cells))

write.csv(summary_tbl, file.path(tbl_dir, "architecture_transfer_summary.csv"), row.names = FALSE)

score_tbl <- seu_inf@meta.data %>%
  filter(!is.na(arch_predicted_label)) %>%
  group_by(arch_predicted_label) %>%
  summarise(
    n_cells = dplyr::n(),
    mean_prediction_score = mean(arch_prediction_score_max, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

write.csv(score_tbl, file.path(tbl_dir, "architecture_transfer_prediction_scores.csv"), row.names = FALSE)

# -----------------------------
# Figures
# -----------------------------
red_use <- pick_plot_reduction(seu_inf)

if (!is.null(red_use)) {
  p1 <- DimPlot(seu_inf, reduction = red_use, group.by = ct_col, raster = TRUE) +
    ggtitle("Original infection labels")
  
  p2 <- DimPlot(seu_inf, reduction = red_use, group.by = "arch_predicted_label", raster = TRUE, na.value = "grey85") +
    ggtitle("Transferred architecture labels")
  
  save_plot(file.path(fig_dir, "Fig35_Architecture_transfer_comparison.pdf"), p1 + p2, w = 16, h = 7)
} else {
  log_msg("No UMAP reduction found for plotting in 13_fib2_reference_mapping.", log_file = log_file)
}

write_legend(
  "Fig35", "Architecture label transfer onto infection dataset",
  hypothesis = "A FIB2-like stromal state should be recoverable inside the broad stromal/perivascular infection labels if the proposed metabolic reservoir is biologically real.",
  methods = "We subset the infection object to stromal, vascular, and trophoblast-adjacent states, then transferred architecture labels from an external placenta reference using Seurat transfer anchors.",
  readout = "Evidence supporting the hypothesis is strongest if broad F/PV states consistently receive FIB2-like transferred labels with nontrivial prediction scores.",
  interpretation_template = "- Concentrated FIB2-like predictions inside F/PV states support a refined stromal reservoir.\n- Diffuse or weak predictions argue for a broader stromal niche rather than a sharply separable FIB2 compartment.",
  outfile = file.path(CFG$dirs$legends, "Fig35_legend.md")
)

log_msg("13_fib2_reference_mapping done.", log_file = log_file)