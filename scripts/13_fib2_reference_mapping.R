#!/usr/bin/env Rscript

source("R/config.R")
source("R/helpers_io.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(qs)
  library(patchwork)
})

log_file <- file.path(CFG$dirs$logs, "13_fib2_reference_mapping.log")
log_msg("13_fib2_reference_mapping starting.", log_file = log_file)

fig_dir <- CFG$dirs$figures
obj_dir <- CFG$dirs$objects
tbl_dir <- CFG$dirs$tables
leg_dir <- CFG$dirs$legends

# thresholds
PRED_SCORE_HIGH <- 0.70
PRED_SCORE_REFINE <- 0.85
CROSS_LINEAGE_SCORE <- 0.95
CROSS_LINEAGE_MARGIN <- 0.10
DIMS <- 1:30

pick_plot_reduction <- function(seu_obj) {
  reds <- names(seu_obj@reductions)
  if ("umap_harmony" %in% reds) return("umap_harmony")
  if ("umap_rna" %in% reds) return("umap_rna")
  if ("X_umap" %in% reds) return("X_umap")
  if ("umap" %in% reds) return("umap")
  NULL
}

safe_qread_or_rds <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  if (grepl("\\.qs$", path)) return(qs::qread(path))
  readRDS(path)
}

prep_for_transfer <- function(obj, npcs = 30, nfeatures = 3000) {
  DefaultAssay(obj) <- "RNA"
  if (!"pca" %in% names(obj@reductions)) {
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- tryCatch(
      FindVariableFeatures(obj, nfeatures = nfeatures, verbose = FALSE),
      error = function(e) {
        log_msg("FindVariableFeatures failed: ", conditionMessage(e), ". Retrying with 2000.", log_file = log_file)
        gc()
        FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
      }
    )
    obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
    obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = npcs, verbose = FALSE)
  }
  obj
}

get_top3_predictions <- function(pred_score_mat) {
  pred_score_mat <- as.matrix(pred_score_mat)
  top3 <- t(apply(pred_score_mat, 2, function(x) {
    ord <- order(x, decreasing = TRUE)[1:min(3, length(x))]
    labs <- rownames(pred_score_mat)[ord]
    vals <- x[ord]

    while (length(labs) < 3) labs <- c(labs, NA_character_)
    while (length(vals) < 3) vals <- c(vals, NA_real_)

    c(
      ref_pred_1 = labs[1], ref_score_1 = vals[1],
      ref_pred_2 = labs[2], ref_score_2 = vals[2],
      ref_pred_3 = labs[3], ref_score_3 = vals[3]
    )
  }))

  top3 <- as.data.frame(top3, stringsAsFactors = FALSE)
  for (cc in c("ref_score_1", "ref_score_2", "ref_score_3")) top3[[cc]] <- as.numeric(top3[[cc]])
  top3
}

map_author_lineage <- function(x) {
  x <- as.character(x)
  out <- rep("other", length(x))
  out[x %in% c("F", "F_p", "F_sm", "PV")] <- "stromal"
  out[x %in% c("Endo_f")] <- "endothelial"
  out[x %in% c("EVT_1", "EVT_2", "iEVT")] <- "EVT"
  out[x %in% c("VCT", "VCT_fusing", "VCT_p", "VCT_CCC")] <- "villous_trophoblast"
  out[x %in% c("HBC", "HBC_p", "PAMM1")] <- "immune_myeloid"
  out
}

map_ref_lineage <- function(x) {
  x <- as.character(x)
  out <- rep("other", length(x))
  out[grepl("FIB|Fibroblast|Maternal\\.fibroblast|Mat\\. ?FIB", x, ignore.case = TRUE)] <- "stromal"
  out[grepl("Endothelial|Endo", x, ignore.case = TRUE)] <- "endothelial"
  out[grepl("^EVT|EVT-progenitor", x, ignore.case = TRUE)] <- "EVT"
  out[grepl("^vCTB|^STB|STB-progenitor|VCT", x, ignore.case = TRUE)] <- "villous_trophoblast"
  out[grepl("Hofbauer|maternal\\.macrophages|MAC|myeloid|PAMM", x, ignore.case = TRUE)] <- "immune_myeloid"
  out
}

is_unknown_label <- function(x) {
  is.na(x) | x %in% c("", "NA", "Unknown", "unknown", "unassigned", "Unassigned")
}

# ------------------
# Load query object
# ------------------
query_candidates <- c(
  file.path(obj_dir, "seu_with_integration_checks.qs"),
  file.path(obj_dir, "seu_with_ea_flt1_proxies.qs"),
  file.path(obj_dir, "seu_with_scores.qs")
)
query_path <- query_candidates[file.exists(query_candidates)][1]
if (is.na(query_path) || length(query_path) == 0) {
  stop("Could not find query object. Checked:\n", paste(query_candidates, collapse = "\n"))
}
log_msg("Using query object: ", query_path, log_file = log_file)

seu_inf <- safe_qread_or_rds(query_path)
DefaultAssay(seu_inf) <- "RNA"
ct_col <- CFG$cols$cell_type
if (!ct_col %in% colnames(seu_inf@meta.data)) stop("Missing query cell type column: ", ct_col)
seu_inf$celltype_author <- as.character(seu_inf@meta.data[[ct_col]])

# ------------------
# Load reference
# ------------------
ref_candidates <- c(
  file.path("data", "02_processed", "slidetags_mapped_to_multiome.rds"),
  file.path("data", "processed", "slidetags_mapped_to_multiome.rds"),
  file.path("outputs", "objects", "slidetags_mapped_to_multiome.qs"),
  file.path("outputs", "objects", "slidetags_mapped_to_multiome.rds"),
  file.path("..", "hPlacenta-architecture", "data", "processed", "slidetags_mapped_to_multiome.rds")
)
ref_path <- ref_candidates[file.exists(ref_candidates)][1]
if (is.na(ref_path) || length(ref_path) == 0) {
  stop("Could not find Slide-tags reference. Checked:\n", paste(ref_candidates, collapse = "\n"))
}
log_msg("Using reference object: ", ref_path, log_file = log_file)

seu_ref <- safe_qread_or_rds(ref_path)
DefaultAssay(seu_ref) <- "RNA"

ref_label_col <- c("cell_type_cluster", "predicted.id")[c("cell_type_cluster", "predicted.id") %in% colnames(seu_ref@meta.data)][1]
if (is.na(ref_label_col) || length(ref_label_col) == 0) stop("Reference label column not found (expected cell_type_cluster or predicted.id)")

# shared genes
shared <- intersect(rownames(seu_inf), rownames(seu_ref))
if (length(shared) < 500) stop("Too few shared genes: ", length(shared))
seu_inf <- subset(seu_inf, features = shared)
seu_ref <- subset(seu_ref, features = shared)

# focus subset for speed
keep_inf <- c("F", "F_p", "F_sm", "PV", "Endo_f", "EVT_1", "EVT_2", "iEVT", "VCT", "VCT_fusing", "VCT_p", "VCT_CCC", "HBC", "HBC_p", "PAMM1")
keep_inf <- intersect(keep_inf, unique(as.character(seu_inf@meta.data[[ct_col]])))
cells_keep_inf <- rownames(seu_inf@meta.data)[as.character(seu_inf@meta.data[[ct_col]]) %in% keep_inf]
if (length(cells_keep_inf) == 0) stop("No query cells matched target lineages.")

seu_inf_sub <- subset(seu_inf, cells = cells_keep_inf)
seu_ref_sub <- seu_ref

seu_inf_sub <- prep_for_transfer(seu_inf_sub, npcs = 30, nfeatures = 3000)
seu_ref_sub <- prep_for_transfer(seu_ref_sub, npcs = 30, nfeatures = 3000)

# transfer
log_msg("Finding transfer anchors...", log_file = log_file)
anchors <- FindTransferAnchors(
  reference = seu_ref_sub,
  query = seu_inf_sub,
  reference.reduction = "pca",
  normalization.method = "LogNormalize",
  dims = DIMS,
  verbose = FALSE
)

log_msg("Running TransferData...", log_file = log_file)
seu_inf_sub <- TransferData(
  anchorset = anchors,
  refdata = list(refined_label = as.character(seu_ref_sub@meta.data[[ref_label_col]])),
  reference = seu_ref_sub,
  query = seu_inf_sub,
  dims = DIMS,
  prediction.assay = TRUE,
  verbose = FALSE
)

pred_label_col <- "predicted.refined_label"
pred_score_col <- "predicted.refined_label.score"
pred_assay_name <- "prediction.score.refined_label"

if (!pred_label_col %in% colnames(seu_inf_sub@meta.data)) stop("Missing prediction label column: ", pred_label_col)
if (!pred_score_col %in% colnames(seu_inf_sub@meta.data)) stop("Missing prediction score column: ", pred_score_col)
if (!pred_assay_name %in% names(seu_inf_sub@assays)) stop("Missing prediction assay: ", pred_assay_name)

pred_mat <- GetAssayData(seu_inf_sub[[pred_assay_name]], slot = "data")
top3 <- get_top3_predictions(pred_mat)
rownames(top3) <- colnames(seu_inf_sub)
seu_inf_sub <- AddMetaData(seu_inf_sub, metadata = top3)

seu_inf_sub$arch_predicted_label <- as.character(seu_inf_sub@meta.data[[pred_label_col]])
seu_inf_sub$arch_prediction_score_max <- as.numeric(seu_inf_sub@meta.data[[pred_score_col]])
seu_inf_sub$prediction_margin_1v2 <- seu_inf_sub$ref_score_1 - seu_inf_sub$ref_score_2

# conservative refinement
seu_inf_sub$author_lineage <- map_author_lineage(seu_inf_sub$celltype_author)
seu_inf_sub$pred_lineage <- map_ref_lineage(seu_inf_sub$ref_pred_1)

same_lineage <- !is_unknown_label(seu_inf_sub$celltype_author) &
  seu_inf_sub$author_lineage == seu_inf_sub$pred_lineage &
  seu_inf_sub$author_lineage != "other"

within_lineage_refine <- same_lineage & seu_inf_sub$ref_score_1 >= PRED_SCORE_REFINE
cross_lineage_high <- !is_unknown_label(seu_inf_sub$celltype_author) &
  seu_inf_sub$author_lineage != seu_inf_sub$pred_lineage &
  seu_inf_sub$ref_score_1 >= CROSS_LINEAGE_SCORE &
  seu_inf_sub$prediction_margin_1v2 >= CROSS_LINEAGE_MARGIN

missing_author_high <- is_unknown_label(seu_inf_sub$celltype_author) &
  seu_inf_sub$ref_score_1 >= PRED_SCORE_HIGH

seu_inf_sub$celltype_refined <- seu_inf_sub$celltype_author
seu_inf_sub$refinement_action <- "kept_author"
seu_inf_sub$celltype_refined[missing_author_high] <- seu_inf_sub$ref_pred_1[missing_author_high]
seu_inf_sub$refinement_action[missing_author_high] <- "filled_missing_with_prediction"
seu_inf_sub$celltype_refined[within_lineage_refine] <- seu_inf_sub$ref_pred_1[within_lineage_refine]
seu_inf_sub$refinement_action[within_lineage_refine] <- "within_lineage_refine"
seu_inf_sub$celltype_refined[cross_lineage_high] <- seu_inf_sub$ref_pred_1[cross_lineage_high]
seu_inf_sub$refinement_action[cross_lineage_high] <- "cross_lineage_override_very_high_confidence"

seu_inf_sub$fib_subtype <- ifelse(
  grepl("FIB1|FIB2|Fibroblast1|Fibroblast2", seu_inf_sub$celltype_refined, ignore.case = TRUE),
  seu_inf_sub$celltype_refined,
  NA_character_
)

seu_inf_sub$immune_identity_refined <- ifelse(
  grepl("Hofbauer|macrophage|myeloid|PAMM", seu_inf_sub$celltype_refined, ignore.case = TRUE),
  seu_inf_sub$celltype_refined,
  NA_character_
)

# map back to full object
new_cols <- c(
  "arch_predicted_label", "arch_prediction_score_max",
  "ref_pred_1", "ref_score_1", "ref_pred_2", "ref_score_2", "ref_pred_3", "ref_score_3",
  "prediction_margin_1v2", "author_lineage", "pred_lineage",
  "celltype_refined", "refinement_action", "fib_subtype", "immune_identity_refined"
)

for (cc in new_cols) seu_inf[[cc]] <- NA
common_cells <- intersect(colnames(seu_inf), colnames(seu_inf_sub))
for (cc in new_cols) seu_inf@meta.data[common_cells, cc] <- seu_inf_sub@meta.data[common_cells, cc]

# save object
qs::qsave(seu_inf, file.path(obj_dir, "seu_with_architecture_transfer.qs"), preset = "high")

# tables
summary_tbl <- seu_inf@meta.data %>%
  filter(!is.na(celltype_refined)) %>%
  count(celltype_author, celltype_refined, refinement_action, name = "n_cells") %>%
  arrange(desc(n_cells))
write.csv(summary_tbl, file.path(tbl_dir, "architecture_transfer_summary.csv"), row.names = FALSE)

score_tbl <- seu_inf@meta.data %>%
  filter(!is.na(arch_predicted_label)) %>%
  group_by(arch_predicted_label) %>%
  summarise(
    n_cells = dplyr::n(),
    mean_prediction_score = mean(as.numeric(arch_prediction_score_max), na.rm = TRUE),
    median_prediction_score = median(as.numeric(arch_prediction_score_max), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))
write.csv(score_tbl, file.path(tbl_dir, "architecture_transfer_prediction_scores.csv"), row.names = FALSE)

thresholds <- c(0.85, 0.90, 0.95)
sens_tbl <- bind_rows(lapply(thresholds, function(t) {
  use_t <- !is.na(seu_inf$ref_score_1) & (seu_inf$ref_score_1 >= t)
  data.frame(
    threshold = t,
    n_refined = sum(use_t, na.rm = TRUE),
    n_cross_lineage = sum(use_t & seu_inf$refinement_action == "cross_lineage_override_very_high_confidence", na.rm = TRUE)
  )
}))
write.csv(sens_tbl, file.path(tbl_dir, "architecture_transfer_threshold_sensitivity.csv"), row.names = FALSE)

# figures
red_use <- pick_plot_reduction(seu_inf)
if (!is.null(red_use)) {
  p1 <- DimPlot(seu_inf, reduction = red_use, group.by = "celltype_author", raster = TRUE) + ggtitle("Original author labels")
  p2 <- DimPlot(seu_inf, reduction = red_use, group.by = "celltype_refined", raster = TRUE, na.value = "grey85") + ggtitle("Refined labels")
  p3 <- FeaturePlot(seu_inf, reduction = red_use, features = "arch_prediction_score_max", raster = TRUE) + ggtitle("Max transfer confidence")

  save_plot(file.path(fig_dir, "Fig35_Architecture_transfer_comparison.pdf"), p1 + p2, w = 16, h = 7)
  save_plot(file.path(fig_dir, "Fig35b_Architecture_transfer_confidence.pdf"), p3, w = 8, h = 6)
} else {
  log_msg("No suitable reduction found for plotting in 13_fib2_reference_mapping.", log_file = log_file)
}

write_legend(
  "Fig35", "Architecture-informed label refinement of infection dataset",
  hypothesis = "Broad stromal/trophoblast/immune labels can be conservatively refined using Slide-tags atlas transfer while keeping author labels as default.",
  methods = "Transfer anchors + TransferData with top-3 prediction capture, conservative lineage-aware refinement thresholds, and back-mapping to full object metadata.",
  readout = "Look for lineage-consistent refinement (e.g., fibro states to FIB1/FIB2) with high confidence and limited cross-lineage overrides.",
  interpretation_template = "- Within-lineage refinements are more reliable than cross-lineage changes.\n- Cross-lineage overrides should be rare and high-confidence.\n- Use threshold sensitivity table to assess robustness.",
  outfile = file.path(leg_dir, "Fig35_legend.md")
)

log_msg("13_fib2_reference_mapping done.", log_file = log_file)
