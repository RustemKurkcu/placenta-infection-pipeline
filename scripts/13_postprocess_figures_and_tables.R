#!/usr/bin/env Rscript
# scripts/13_mapping/13_postprocess_figures_and_tables.R
#
# Postprocessing for Script 13:
# - canonicalize labels for plotting
# - produce additional tables: confusion matrices, per-author summaries, top markers
# - produce additional figures: confusion heatmap, violin plot of confidence, refinement barplot,
#   DotPlot of marker genes, fibro UMAP/marker plots
# - save enhanced figure legends and a Methods/Script-13 md
#
# Run after 13_fib2_reference_mapping.R or with seu_inf loaded.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(qs)
  library(scales)
  library(tidyr)
  library(stringr)
  library(patchwork)
  # forcats used in plotting; reference explicitly when needed
  library(forcats)
})

# ---- Config: ensure these variables are set or override here ----
# If you run interactively and CFG is available, use it; otherwise set paths manually
if (exists("CFG") && is.list(CFG) && !is.null(CFG$dirs)) {
  fig_dir <- CFG$dirs$figures
  tbl_dir <- CFG$dirs$tables
  obj_dir <- CFG$dirs$objects
  leg_dir <- CFG$dirs$legends
} else {
  # adjust these to your repo layout if needed:
  obj_dir <- "outputs/objects"
  tbl_dir <- "outputs/tables"
  fig_dir <- "outputs/figures"
  leg_dir <- "outputs/legends"
}

# Make sure folders exist
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(leg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("scripts/13_mapping", recursive = TRUE, showWarnings = FALSE) # for methods file

safe_msg <- function(...) message(paste0("[13_postprocess] ", paste(..., collapse = " ")))

# ---- Load object if not in environment ----
if (!exists("seu_inf")) {
  # try to read the object we saved earlier
  candidate <- file.path(obj_dir, "seu_with_architecture_transfer.qs")
  if (file.exists(candidate)) {
    safe_msg("Loading Seurat object from: ", candidate)
    seu_inf <- qs::qread(candidate)
  } else {
    stop("seu_inf is not in environment and no saved object found at ", candidate,
         "\nEither load your infection Seurat object as 'seu_inf' or set obj_dir and rerun.")
  }
}
DefaultAssay(seu_inf) <- "RNA"

# ---- Helper functions ----
# Robust canonicalize_label: handles NAs safely and maps exact matches (case-insensitive)
canonicalize_label <- function(x) {
  # Convert to character and preserve NA
  x0 <- as.character(x)
  x1 <- ifelse(is.na(x0), NA_character_, x0)
  
  # basic normalization
  x1 <- gsub("\\.", " ", x1)   # dots -> space
  x1 <- gsub("_", " ", x1)     # underscores -> space
  x1 <- stringr::str_squish(x1) # collapse repeated spaces & trim
  
  # mapping table (extend or change as needed)
  map <- c(
    "hofbauer cells" = "Hofbauer cells",
    "hofbauer.cells" = "Hofbauer cells",
    "hofbauer" = "Hofbauer cells",
    "fibroblast1" = "FIB1",
    "fibroblast2" = "FIB2",
    "fibroblast 1" = "FIB1",
    "fibroblast 2" = "FIB2",
    "fibroblast" = "FIB",
    "maternal fibroblast" = "Maternal.fibroblast",
    "mat fib" = "Mat.FIB",
    "mat.fib" = "Mat.FIB",
    "evt progenitor" = "EVT-progenitor",
    "vctb" = "vCTB",
    "vct" = "VCT",
    "stb progenitor" = "STB-progenitor"
  )
  
  # prepare keys & values
  map_keys <- names(map)
  map_vals <- unname(map)
  
  # lowercased vector for matching
  x1_l <- tolower(x1)
  
  non_na_inds <- which(!is.na(x1_l))
  if (length(non_na_inds) > 0) {
    # match exact lowercased strings
    matched <- match(x1_l[non_na_inds], map_keys)
    valid <- which(!is.na(matched))
    if (length(valid) > 0) {
      assign_inds <- non_na_inds[valid]
      x1[assign_inds] <- map_vals[ matched[valid] ]
    }
  }
  
  x1
}

# safe save plot helper using ggsave or grid/pdf
save_plot_safe <- function(filename, plot_obj, width=8, height=6) {
  try({
    if (inherits(plot_obj, "ggplot") || inherits(plot_obj, "patchwork")) {
      ggsave(filename, plot = plot_obj, width = width, height = height, units = "in", device = "pdf")
    } else {
      # assume plot_obj is a grid object
      pdf(filename, width = width, height = height)
      print(plot_obj)
      dev.off()
    }
  }, silent = FALSE)
}

# ---- Canonicalize labels in metadata for plots and tables ----
safe_msg("Canonicalizing label columns for consistent plotting/tables...")
label_cols <- c("arch_predicted_label","pred1","pred2","pred3",
                "celltype_refined","celltype_author","refinement_action",
                "fib_subtype","immune_identity_refined")
label_cols <- intersect(label_cols, colnames(seu_inf@meta.data))
for (cc in label_cols) {
  safe_msg(" - canonicalizing column:", cc)
  seu_inf@meta.data[[cc]] <- canonicalize_label(seu_inf@meta.data[[cc]])
}
safe_msg("Canonicalization done.")

# ---- 1) Confusion matrix: author -> refined (counts and pct by author) ----
if (!("celltype_author" %in% colnames(seu_inf@meta.data)) || !("celltype_refined" %in% colnames(seu_inf@meta.data))) {
  stop("Required columns 'celltype_author' or 'celltype_refined' are missing from metadata.")
}

safe_msg("Building confusion matrix (author -> refined)...")
conf_counts <- table(seu_inf@meta.data$celltype_author, seu_inf@meta.data$celltype_refined, dnn = c("author","refined"))
conf_counts_df <- as.data.frame.matrix(conf_counts)
write.csv(conf_counts_df, file.path(tbl_dir, "confusion_author_to_refined_counts.csv"), row.names = TRUE)
# percents by author (row-wise)
conf_pct_by_author <- prop.table(conf_counts, margin=1)
conf_pct_df <- as.data.frame.matrix(round(conf_pct_by_author,3))
write.csv(conf_pct_df, file.path(tbl_dir, "confusion_author_to_refined_pct_by_author.csv"), row.names = TRUE)
safe_msg("Wrote confusion tables.")

# heatmap (normalized by author)
safe_msg("Plotting confusion heatmap (rows=author, cols=refined)...")
try({
  pngfile <- file.path(fig_dir, "Fig36_confusion_heatmap.png")
  pheatmap::pheatmap(as.matrix(conf_pct_by_author),
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     main = "Author -> Refined (proportion by author)",
                     fontsize = 9,
                     filename = pngfile,
                     width = 9, height = 7)
  safe_msg("Wrote heatmap to:", pngfile)
}, silent = FALSE)

# ---- 2) Per-author summary table ----
safe_msg("Building per-author summary table...")
per_author <- seu_inf@meta.data %>%
  group_by(celltype_author) %>%
  summarise(
    n_cells = n(),
    n_refined = sum(!is.na(celltype_refined) & celltype_refined != celltype_author),
    prop_refined = ifelse(n_cells>0, n_refined / n_cells, NA_real_),
    mean_ref_score = mean(as.numeric(arch_prediction_score_max), na.rm = TRUE),
    median_ref_score = median(as.numeric(arch_prediction_score_max), na.rm = TRUE),
    n_cross_lineage_overrides = sum(refinement_action == "cross_lineage_override_very_high_confidence", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))
write.csv(per_author, file.path(tbl_dir, "per_author_summary.csv"), row.names = FALSE)
safe_msg("Wrote per_author_summary.csv")

# ---- 3) Violin plot: arch_prediction_score_max by celltype_author ----
safe_msg("Creating violin plot of arch_prediction_score_max by author...")
if ("arch_prediction_score_max" %in% colnames(seu_inf@meta.data)) {
  df_v <- seu_inf@meta.data %>%
    mutate(arch_prediction_score_max = as.numeric(arch_prediction_score_max)) %>%
    filter(!is.na(celltype_author))
  p_violin <- ggplot(df_v, aes(x = forcats::fct_reorder(celltype_author, arch_prediction_score_max, .fun = median, .desc = TRUE), y = arch_prediction_score_max)) +
    geom_violin(fill="#2b83ba", alpha=0.6) +
    geom_boxplot(width=0.12, outlier.size=0.5) +
    coord_flip() +
    theme_bw() + labs(x = "Author cell type", y = "Max architecture prediction score", title = "Transfer confidence by author cell type") +
    theme(axis.text.y = element_text(size=8))
  ggsave(file.path(fig_dir, "Fig37_prediction_score_by_author.pdf"), p_violin, width = 10, height = 8)
  safe_msg("Wrote Fig37_prediction_score_by_author.pdf")
} else {
  safe_msg("arch_prediction_score_max not present; skipping violin plot.")
}

# ---- 4) Refinement action bar plot by author ----
safe_msg("Creating refinement action barplot by author...")
df_bar <- seu_inf@meta.data %>%
  mutate(celltype_author = as.character(celltype_author),
         refinement_action = as.character(refinement_action)) %>%
  group_by(celltype_author, refinement_action) %>%
  summarize(n = n(), .groups="drop") %>%
  group_by(celltype_author) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
p_bar <- ggplot(df_bar, aes(x = forcats::fct_reorder(celltype_author, prop, .fun = sum), y = prop, fill=refinement_action)) +
  geom_col() + coord_flip() + theme_minimal() +
  labs(x = "Author cell type", y = "Proportion", fill = "Refinement action", title = "Refinement actions per author cell type")
ggsave(file.path(fig_dir, "Fig38_refinement_action_by_author.pdf"), p_bar, width = 10, height = 8)
safe_msg("Wrote Fig38_refinement_action_by_author.pdf")

# ---- 5) Top-3 table already produced by main script; ensure it's present, otherwise produce minimal version ----
if (!file.exists(file.path(tbl_dir, "architecture_transfer_top3_predictions.csv"))) {
  safe_msg("Top3 CSV not found; creating from metadata.")
  top3_tbl <- seu_inf@meta.data %>%
    mutate(cell = rownames(seu_inf@meta.data)) %>%
    filter(!is.na(pred1)) %>%
    transmute(cell, celltype_author, pred1, score1 = as.numeric(score1),
              pred2, score2 = as.numeric(score2), pred3, score3 = as.numeric(score3),
              prediction_margin_1v2 = as.numeric(prediction_margin_1v2),
              celltype_refined, refinement_action)
  write.csv(top3_tbl, file.path(tbl_dir, "architecture_transfer_top3_predictions.csv"), row.names = FALSE)
  safe_msg("Wrote architecture_transfer_top3_predictions.csv")
} else {
  safe_msg("Top3 CSV already exists; skipped regeneration.")
}

# ---- 6) Threshold sensitivity (again, but saved) ----
safe_msg("Recreating threshold sensitivity table...")
thresholds <- c(0.85, 0.90, 0.95)
sens_tbl <- lapply(thresholds, function(t) {
  use_t <- !is.na(seu_inf@meta.data$score1) & (as.numeric(seu_inf@meta.data$score1) >= t)
  tibble::tibble(threshold = t,
                 n_overrides = sum(use_t & (seu_inf@meta.data$pred1 != seu_inf@meta.data$celltype_author), na.rm = TRUE),
                 n_high_conf = sum(use_t, na.rm = TRUE))
})
sens_tbl <- dplyr::bind_rows(sens_tbl)
write.csv(sens_tbl, file.path(tbl_dir, "architecture_transfer_threshold_sensitivity.csv"), row.names = FALSE)
safe_msg("Wrote architecture_transfer_threshold_sensitivity.csv")

# ---- 7) Top markers for fibroblast subtypes ----
safe_msg("Finding top markers for fibroblast subtypes (if labels present)...")
# Guard: ensure celltype_refined exists and is a valid identity
if ("celltype_refined" %in% colnames(seu_inf@meta.data) &&
    any(grepl("^FIB", toupper(as.character(seu_inf@meta.data$celltype_refined)), perl = TRUE))) {
  
  # set idents to celltype_refined safely
  seu_inf <- Seurat::SetIdent(seu_inf, value = seu_inf@meta.data$celltype_refined)
  refined_labels <- unique(as.character(seu_inf@meta.data$celltype_refined))
  
  if (all(c("FIB2","FIB1") %in% toupper(refined_labels))) {
    lab1 <- refined_labels[toupper(refined_labels) == "FIB2"][1]
    lab2 <- refined_labels[toupper(refined_labels) == "FIB1"][1]
    if (!is.na(lab1) && !is.na(lab2)) {
      try({
        markers_fib2_vs_fib1 <- FindMarkers(seu_inf, ident.1 = lab1, ident.2 = lab2, min.pct = 0.1, logfc.threshold = 0.25)
        markers_fib2_vs_fib1 <- markers_fib2_vs_fib1 %>% tibble::rownames_to_column("gene") %>% arrange(p_val_adj, desc(avg_log2FC))
        write.csv(markers_fib2_vs_fib1, file.path(tbl_dir, "FIB2_vs_FIB1_markers.csv"), row.names = FALSE)
        safe_msg("Wrote FIB2_vs_FIB1_markers.csv")
      }, silent = FALSE)
    }
  } else {
    try({
      markers_all <- FindAllMarkers(seu_inf, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
      write.csv(markers_all, file.path(tbl_dir, "markers_by_refined_label.csv"), row.names = FALSE)
      safe_msg("Wrote markers_by_refined_label.csv")
    }, silent = FALSE)
  }
} else {
  safe_msg("No FIB1/FIB2 labels detected; skipping fibro marker finding.")
}

# ---- 8) DotPlot for canonical FIB2 markers and selected markers across refined labels ----
safe_msg("Generating DotPlot for canonical FIB2/stromal markers and trophoblast/immune markers...")
markers_candidate <- c("PDGFRA","CXCL14","DCN","SPON1","DLK1","COL6A3","AREG","HGF","VEGFA",
                       "IL1B","TNF","NOS2","CD163","MRC1","ARG1",
                       "GATA3","KRT7","TP63")
present_genes <- intersect(markers_candidate, rownames(seu_inf))
if (length(present_genes) >= 3 && "celltype_refined" %in% colnames(seu_inf@meta.data)) {
  try({
    dp <- DotPlot(seu_inf, features = present_genes, group.by = "celltype_refined") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle("Selected marker DotPlot across refined labels")
    ggsave(file.path(fig_dir, "Fig39_marker_DotPlot_by_refined.pdf"), dp, width = 12, height = 6)
    safe_msg("Wrote Fig39_marker_DotPlot_by_refined.pdf")
  }, silent = FALSE)
} else {
  safe_msg("Fewer than 3 candidate genes present or celltype_refined missing; skipping DotPlot.")
}

# ---- 9) Subset UMAP for fibro/stromal with markers featureplots (if UMAP exists) ----
# Choose reduction preference
red_use <- NULL
for (cand in c("umap_harmony","umap_rna","X_umap","umap","umap_harm")) {
  if (cand %in% names(seu_inf@reductions)) { red_use <- cand; break }
}
safe_msg("Using reduction:", red_use)
if (!is.null(red_use) && red_use %in% names(seu_inf@reductions)) {
  fibro_mask <- FALSE
  if ("celltype_refined" %in% colnames(seu_inf@meta.data)) {
    fibro_mask <- grepl("^FIB|Fibro|^F\\b|PV|F_p|F_sm", seu_inf@meta.data$celltype_refined, ignore.case = TRUE)
  }
  if (!any(fibro_mask) && "celltype_author" %in% colnames(seu_inf@meta.data)) {
    fibro_mask <- grepl("^FIB|Fibro|^F\\b|PV|F_p|F_sm", seu_inf@meta.data$celltype_author, ignore.case = TRUE)
  }
  fibro_cells <- colnames(seu_inf)[which(fibro_mask)]
  if (length(fibro_cells) > 10) {
    seu_fib <- subset(seu_inf, cells = fibro_cells)
    fibro_features <- c("PDGFRA","DCN","CXCL14","SPON1")
    fibro_features_present <- intersect(fibro_features, rownames(seu_fib))
    if (length(fibro_features_present) > 0) {
      try({
        fp <- FeaturePlot(seu_fib, features = fibro_features_present, reduction = red_use, pt.size = 0.5, ncol = 2)
        save_plot_safe(file.path(fig_dir, "Fig40_fibro_featureplots.pdf"), fp, width = 10, height = 8)
        safe_msg("Wrote Fig40_fibro_featureplots.pdf")
      }, silent = FALSE)
    } else {
      safe_msg("No fibro canonical markers found in object for FeaturePlot.")
    }
  } else {
    safe_msg("Too few fibro-like cells for fibro-specific UMAP/plots.")
  }
} else {
  safe_msg("No suitable reduction found for fibro UMAP plotting.")
}

# ---- 10) Write enhanced figure legend (markdown) and methods file ----
safe_msg("Writing enhanced figure legend and Script 13 methods file...")
fig35_legend <- c(
  "# Fig35 — Architecture informed label refinement",
  "",
  "**Panel A — Author labels.** UMAP projection of the infection dataset colored by the original author-supplied cell types (`celltype_author`).",
  "",
  "**Panel B — Refined labels after architecture transfer.** The same UMAP colored by `celltype_refined` after conservative architecture-based label transfer from the Slide-tags mapped multiome reference. Refinement rules applied: author labels are retained by default; within-lineage refinements occur at `score >= 0.85`; cross-lineage overrides allowed only at `score >= 0.95` and `margin >= 0.10`.",
  "",
  "**Panel C — Transfer confidence.** Continuous color scale displays maximum transfer confidence (`arch_prediction_score_max`).",
  "",
  "**Interpretation notes:**",
  "- Concentration of FIB2 predictions within stromal author labels supports the fibroblast reservoir hypothesis.",
  "- Consult `architecture_transfer_top3_predictions.csv` for per-cell audit."
)
writeLines(fig35_legend, con = file.path(leg_dir, "Fig35_legend.md"))

methods_text <- c(
  "# Script 13 methods — Architecture-based label transfer and conservative refinement",
  "",
  "**Objective.** Refine broad author labels in the infection dataset using a high-resolution placenta Slide-tags reference, focusing on robust identification of fibroblast subtypes (FIB1/FIB2) and improved immune/myeloid annotation.",
  "",
  "**Approach.** We compute PCA on query and reference (only if needed), build transfer anchors using Seurat `FindTransferAnchors()` (dims=1:30), then `TransferData()` (k.weight=50, weight.reduction='pcaproject') to obtain predicted labels and per-label score assays. Top-3 predicted labels are extracted per cell and margins computed.",
  "",
  "**Conservative refinement rules.**",
  "- Keep author labels by default.",
  "- Fill missing author labels if `score1 >= 0.70`.",
  "- Within-lineage refinement if `author_lineage == pred_lineage` and `score1 >= 0.85`.",
  "- Cross-lineage override only if `score1 >= 0.95` and `(score1 - score2) >= 0.10`.",
  "",
  "**Outputs.** Key CSVs: `architecture_transfer_summary.csv`, `architecture_transfer_top3_predictions.csv`, `architecture_transfer_prediction_scores.csv`, `architecture_transfer_threshold_sensitivity.csv`, `refinement_action_by_author.csv`.",
  "",
  "**Troubleshooting.** If you hit memory errors: reduce DIMS, reduce k.weight, or run on a larger node."
)
writeLines(methods_text, con = file.path("scripts", "13_mapping", "13_methods.md"))

safe_msg("Legend and methods files written to:", leg_dir, "and scripts/13_mapping respectively.")

# ---- final save of object (with canonicalized labels) ----
safe_msg("Saving canonicalized Seurat object to outputs...")
out_qs <- file.path(obj_dir, "seu_with_architecture_transfer_canonicalized.qs")
qs::qsave(seu_inf, out_qs, preset = "high")
safe_msg("Saved object:", out_qs)

safe_msg("Post-processing complete. Tables in:", tbl_dir, "Figures in:", fig_dir, "Legends in:", leg_dir)