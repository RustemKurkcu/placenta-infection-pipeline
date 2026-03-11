source("R/config.R")
source("R/helpers_io.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(qs)
  library(edgeR)
})

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("10_ea_flt1_pseudobulk starting.", log_file = log_file)

seu <- qs::qread(file.path(CFG$dirs$objects, "seu_with_scores.qs"))
DefaultAssay(seu) <- "RNA"

# --------------------------------------------------
# WHY:
# We already have cell types annotated. Do not redo annotation.
# We aggregate to pseudobulk because donor/sample is the inferential unit.
# --------------------------------------------------

# sample_id
if (!"sample_id" %in% colnames(seu@meta.data)) {
  seu$sample_id <- paste(
    seu[[CFG$cols$donor_id]][,1],
    seu[[CFG$cols$infection]][,1],
    seu[[CFG$cols$hpi]][,1],
    seu[[CFG$cols$stage]][,1],
    sep = "|"
  )
}

# Combined condition label
if (!"condition" %in% colnames(seu@meta.data)) {
  seu$condition <- paste0(seu[[CFG$cols$infection]][,1], "_", seu[[CFG$cols$hpi]][,1])
}

# A unique pseudobulk ID: sample + cell type
seu$pb_id <- paste(seu$sample_id, seu[[CFG$cols$cell_type]][,1], sep = "||")

# --------------------------------------------------
# Curated gene sets
# NOTE:
# FLT1 here is TOTAL FLT1 transcript signal, not sFlt1-e15a specifically.
# Isoform-specific analysis needs exon/junction-aware data or targeted assays.
# --------------------------------------------------

ea_release_genes <- c("PLD1", "PLD2", "GDPD1", "GDPD5", "FAAH", "NAPEPLD")
ea_recycle_genes <- c("ETNK1", "ETNK2", "PCYT2", "SELENOI", "CEPT1")

flt1_axis_genes  <- c("FLT1", "PGF", "KDR", "VEGFA", "ENG", "HIF1A", "EPAS1")
inflam_genes     <- c("IL1B", "CCL20", "CCL4", "CXCL8", "PTGS2", "LGALS3", "AREG", "HMOX1")

present <- function(x) intersect(x, rownames(seu))
ea_release_genes <- present(ea_release_genes)
ea_recycle_genes <- present(ea_recycle_genes)
flt1_axis_genes  <- present(flt1_axis_genes)
inflam_genes     <- present(inflam_genes)

# --------------------------------------------------
# Pseudobulk counts
# --------------------------------------------------
pb_counts <- AggregateExpression(
  object = seu,
  assays = "RNA",
  group.by = "pb_id",
  slot = "counts",
  return.seurat = FALSE
)$RNA

pb_meta <- seu@meta.data %>%
  distinct(pb_id,
           sample_id,
           donor_id = .data[[CFG$cols$donor_id]],
           infection = .data[[CFG$cols$infection]],
           hpi = .data[[CFG$cols$hpi]],
           stage = .data[[CFG$cols$stage]],
           condition,
           cell_type = .data[[CFG$cols$cell_type]]) %>%
  filter(pb_id %in% colnames(pb_counts))

# Keep only pseudobulks with enough cells
pb_ncells <- seu@meta.data %>%
  count(pb_id, name = "n_cells_pb")

pb_meta <- pb_meta %>%
  left_join(pb_ncells, by = "pb_id") %>%
  filter(n_cells_pb >= 30)

pb_counts <- pb_counts[, pb_meta$pb_id, drop = FALSE]

# --------------------------------------------------
# Normalize pseudobulk
# --------------------------------------------------
dge <- DGEList(counts = pb_counts)
dge <- calcNormFactors(dge)
logcpm <- cpm(dge, log = TRUE, prior.count = 2)

score_mean <- function(mat, genes) {
  genes <- intersect(genes, rownames(mat))
  if (length(genes) == 0) return(rep(NA_real_, ncol(mat)))
  colMeans(mat[genes, , drop = FALSE], na.rm = TRUE)
}

pb_scores <- tibble(
  pb_id = colnames(logcpm),
  EA_release = score_mean(logcpm, ea_release_genes),
  EA_recycle = score_mean(logcpm, ea_recycle_genes),
  FLT1_axis  = score_mean(logcpm, flt1_axis_genes),
  Inflamm    = score_mean(logcpm, inflam_genes),
  FLT1_expr  = if ("FLT1" %in% rownames(logcpm)) logcpm["FLT1", ] else NA_real_,
  PGF_expr   = if ("PGF"  %in% rownames(logcpm)) logcpm["PGF", ]  else NA_real_,
  IL1B_expr  = if ("IL1B" %in% rownames(logcpm)) logcpm["IL1B", ] else NA_real_,
  HMOX1_expr = if ("HMOX1" %in% rownames(logcpm)) logcpm["HMOX1", ] else NA_real_
) %>%
  mutate(
    EA_proxy = EA_release - EA_recycle,
    FLT1_PGF_diff = FLT1_expr - PGF_expr
  ) %>%
  left_join(pb_meta, by = "pb_id")

ensure_dir(CFG$dirs$tables)
write.csv(pb_scores, file.path(CFG$dirs$tables, "pseudobulk_module_scores.csv"), row.names = FALSE)

# --------------------------------------------------
# Donor-aware fixed-effect model
# WHY:
# donor is the replicate unit; infection/time are compared within donor structure
# --------------------------------------------------
fit_feature_by_celltype <- function(df, response) {
  out <- lapply(split(df, df$cell_type), function(dd) {
    # require enough samples
    if (nrow(dd) < 6 || length(unique(dd$infection)) < 2) return(NULL)
    
    dd$donor_id  <- factor(dd$donor_id)
    dd$infection <- factor(dd$infection)
    dd$hpi       <- factor(dd$hpi)
    
    formula_txt <- paste(response, "~ donor_id + infection * hpi")
    fit <- lm(as.formula(formula_txt), data = dd)
    aov_tab <- anova(fit)
    
    data.frame(
      cell_type = unique(dd$cell_type),
      term = rownames(aov_tab),
      df = aov_tab$Df,
      F_value = aov_tab$`F value`,
      p_value = aov_tab$`Pr(>F)`
    )
  })
  
  bind_rows(out)
}

ea_stats    <- fit_feature_by_celltype(pb_scores, "EA_proxy")
flt1_stats  <- fit_feature_by_celltype(pb_scores, "FLT1_PGF_diff")
inflam_stats <- fit_feature_by_celltype(pb_scores, "Inflamm")

write.csv(ea_stats,   file.path(CFG$dirs$tables, "EA_proxy_stats_by_celltype.csv"), row.names = FALSE)
write.csv(flt1_stats, file.path(CFG$dirs$tables, "FLT1_PGF_stats_by_celltype.csv"), row.names = FALSE)
write.csv(inflam_stats, file.path(CFG$dirs$tables, "Inflamm_stats_by_celltype.csv"), row.names = FALSE)

# --------------------------------------------------
# Figures
# --------------------------------------------------
plot_heat <- function(df, value_col, title_txt, outname) {
  agg <- df %>%
    group_by(cell_type, condition) %>%
    summarise(value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")
  
  p <- ggplot(agg, aes(x = condition, y = cell_type, fill = value)) +
    geom_tile() +
    theme_bw() +
    labs(title = title_txt, x = NULL, y = NULL)
  
  save_plot(file.path(CFG$dirs$figures, outname), p, w = 10, h = 7)
}

plot_spaghetti <- function(df, value_col, title_txt, outname, keep_celltypes = NULL) {
  dd <- df
  if (!is.null(keep_celltypes)) dd <- dd %>% filter(cell_type %in% keep_celltypes)
  
  p <- ggplot(dd, aes(x = condition, y = .data[[value_col]], group = donor_id, color = donor_id)) +
    geom_line(alpha = 0.6) +
    geom_point(size = 2) +
    facet_wrap(~cell_type, scales = "free_y") +
    theme_bw() +
    labs(title = title_txt, x = NULL, y = value_col)
  
  save_plot(file.path(CFG$dirs$figures, outname), p, w = 14, h = 9)
}

plot_heat(pb_scores, "EA_proxy", "EA proxy by cell type and condition", "Fig21_EA_proxy_heatmap.pdf")
plot_heat(pb_scores, "FLT1_PGF_diff", "FLT1-PGF pseudobulk difference by cell type and condition", "Fig22_FLT1_PGF_heatmap.pdf")
plot_spaghetti(pb_scores, "EA_proxy", "EA proxy by donor across conditions", "Fig23_EA_proxy_spaghetti.pdf",
               keep_celltypes = c("F", "F_p", "F_sm", "PV", "VCT", "VCT_fusing", "Endo_f"))
plot_spaghetti(pb_scores, "FLT1_PGF_diff", "FLT1-PGF by donor across conditions", "Fig24_FLT1_PGF_spaghetti.pdf",
               keep_celltypes = c("VCT", "VCT_fusing", "Endo_f", "PV"))

# UMAP overlays on existing reductions
score_to_cells <- pb_scores %>%
  select(pb_id, EA_proxy, FLT1_expr, PGF_expr, IL1B_expr, HMOX1_expr)

cell_md <- seu@meta.data %>%
  select(pb_id) %>%
  left_join(score_to_cells, by = "pb_id")

seu$EA_proxy_pb  <- cell_md$EA_proxy
seu$FLT1_pb      <- cell_md$FLT1_expr
seu$PGF_pb       <- cell_md$PGF_expr
seu$IL1B_pb      <- cell_md$IL1B_expr
seu$HMOX1_pb     <- cell_md$HMOX1_expr

red_use <- if ("umap_harmony" %in% names(seu@reductions)) "umap_harmony" else "X_umap"

for (feat in c("EA_proxy_pb", "FLT1_pb", "PGF_pb", "IL1B_pb", "HMOX1_pb")) {
  p <- FeaturePlot(seu, features = feat, reduction = red_use, split.by = "condition", raster = TRUE)
  save_plot(file.path(CFG$dirs$figures, paste0("Fig_", feat, "_split_condition.pdf")), p, w = 16, h = 10)
}

qs::qsave(pb_scores, file.path(CFG$dirs$objects, "pseudobulk_module_scores.qs"))
log_msg("10_ea_flt1_pseudobulk done.", log_file = log_file)