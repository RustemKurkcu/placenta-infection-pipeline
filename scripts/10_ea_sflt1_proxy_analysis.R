source("R/config.R")
source("R/helpers_io.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(qs)
  library(Matrix)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(patchwork)
  library(tidyr)
})

# ============================================================
# 10_ea_sflt1_proxy_analysis.R
#
# HYPOTHESIS
# ----------
# Specific placental cell types, especially stromal/trophoblast/vascular
# compartments, should show condition-dependent shifts in host programs that
# increase ethanolamine availability and/or vascular stress.
#
# IMPORTANT INTERPRETATION
# ------------------------
# This object contains HOST transcriptomes, not bacterial transcripts.
# Therefore these are host-side PROXIES, not direct measurements of bacterial
# EA catabolism, MegL activity, H2S, NH3, or the placental e15a isoform.
# ============================================================

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("10_ea_sflt1_proxy_analysis starting.", log_file = log_file)

fig_dir <- CFG$dirs$figures
tbl_dir <- CFG$dirs$tables
obj_dir <- CFG$dirs$objects

seu <- qs::qread(file.path(obj_dir, "seu_with_scores.qs"))
DefaultAssay(seu) <- "RNA"

if (!"condition" %in% colnames(seu@meta.data) && all(c(CFG$cols$infection, CFG$cols$hpi) %in% colnames(seu@meta.data))) {
  seu$condition <- paste0(seu[[CFG$cols$infection]][,1], "_", seu[[CFG$cols$hpi]][,1])
}

present <- function(x) intersect(x, rownames(seu))

EA_LIBERATION_GENES <- present(c("PLD1", "PLD2", "GDPD1", "GDPD5", "FAAH", "NAPEPLD"))
EA_SINK_GENES       <- present(c("ETNK1", "ETNK2", "PCYT2", "SELENOI", "CEPT1"))
ANGIO_STRESS_GENES  <- present(c("FLT1", "PGF", "ENG", "KDR", "VEGFA", "HIF1A", "EPAS1", "HMOX1", "IL1B", "CXCL8", "CCL20", "PTGS2", "ICAM1", "VCAM1"))

score_mean <- function(seu_obj, genes, slot = "data") {
  genes <- intersect(genes, rownames(seu_obj))
  if (length(genes) == 0) return(rep(NA_real_, ncol(seu_obj)))
  mat <- GetAssayData(seu_obj, assay = DefaultAssay(seu_obj), slot = slot)
  Matrix::colMeans(mat[genes, , drop = FALSE])
}

safe_scale <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

seu$EA_release_score  <- score_mean(seu, EA_LIBERATION_GENES)
seu$EA_recycle_score  <- score_mean(seu, EA_SINK_GENES)
seu$EA_proxy          <- safe_scale(seu$EA_release_score) - safe_scale(seu$EA_recycle_score)
seu$Angio_stress_score <- score_mean(seu, ANGIO_STRESS_GENES)

flt1_expr <- if ("FLT1" %in% rownames(seu)) as.numeric(GetAssayData(seu, slot = "data")["FLT1", ]) else rep(NA_real_, ncol(seu))
pgf_expr  <- if ("PGF"  %in% rownames(seu)) as.numeric(GetAssayData(seu, slot = "data")["PGF", ])  else rep(NA_real_, ncol(seu))
seu$FLT1_expr <- flt1_expr
seu$PGF_expr  <- pgf_expr
seu$FLT1_minus_PGF <- safe_scale(flt1_expr) - safe_scale(pgf_expr)

if (!"pca" %in% names(seu@reductions)) {
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = 4000, verbose = FALSE)
  seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
  seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 50, verbose = FALSE)
}
if (!"umap_rna" %in% names(seu@reductions)) {
  seu <- RunUMAP(seu, dims = 1:40, reduction = "pca", reduction.name = "umap_rna", verbose = FALSE)
}
if (!"tsne_rna" %in% names(seu@reductions) && requireNamespace("Rtsne", quietly = TRUE)) {
  seu <- RunTSNE(seu, dims = 1:40, reduction = "pca", reduction.name = "tsne_rna", check_duplicates = FALSE, verbose = FALSE)
}

plot_score <- function(feature, title_prefix) {
  p1 <- FeaturePlot(seu, features = feature, reduction = "umap_rna", cols = c("white", "red"), min.cutoff = "q05", max.cutoff = "q95") +
    ggtitle(paste0(title_prefix, " | UMAP"))
  plots <- list(p1)
  if ("tsne_rna" %in% names(seu@reductions)) {
    p2 <- FeaturePlot(seu, features = feature, reduction = "tsne_rna", cols = c("white", "red"), min.cutoff = "q05", max.cutoff = "q95") +
      ggtitle(paste0(title_prefix, " | tSNE"))
    plots <- list(p1, p2)
  }
  wrap_plots(plots)
}

save_plot(file.path(fig_dir, "Fig21_EA_proxy_UMAP_tSNE.pdf"), plot_score("EA_proxy", "EA proxy"), w = 16, h = 7)
save_plot(file.path(fig_dir, "Fig22_AngioStress_UMAP_tSNE.pdf"), plot_score("Angio_stress_score", "Angiogenic stress score"), w = 16, h = 7)
save_plot(file.path(fig_dir, "Fig23_FLT1minusPGF_UMAP_tSNE.pdf"), plot_score("FLT1_minus_PGF", "FLT1 minus PGF proxy"), w = 16, h = 7)

split_plot <- function(feature, title_txt) {
  FeaturePlot(seu, features = feature, split.by = "condition", reduction = "umap_rna", cols = c("white", "red"), min.cutoff = "q05", max.cutoff = "q95") +
    ggtitle(title_txt)
}

save_plot(file.path(fig_dir, "Fig24_EA_proxy_split_condition.pdf"), split_plot("EA_proxy", "EA proxy by condition"), w = 18, h = 10)
save_plot(file.path(fig_dir, "Fig25_FLT1minusPGF_split_condition.pdf"), split_plot("FLT1_minus_PGF", "FLT1 minus PGF proxy by condition"), w = 18, h = 10)

key_genes <- intersect(c("PLD1", "PLD2", "GDPD1", "GDPD5", "FAAH", "NAPEPLD", "ETNK1", "ETNK2", "PCYT2", "SELENOI", "CEPT1", "FLT1", "PGF", "HMOX1", "ENG", "IL1B", "CXCL8", "CCL20"), rownames(seu))
if (length(key_genes) > 0) {
  p_dot <- DotPlot(seu, features = key_genes, group.by = CFG$cols$cell_type) + RotatedAxis() + ggtitle("EA / FLT1-axis genes by cell type")
  save_plot(file.path(fig_dir, "Fig26_Dotplot_EA_FLT1_genes_by_celltype.pdf"), p_dot, w = 16, h = 7)
}

md <- seu@meta.data %>%
  mutate(
    donor_id = .data[[CFG$cols$donor_id]],
    infection = .data[[CFG$cols$infection]],
    hpi = .data[[CFG$cols$hpi]],
    cell_type = .data[[CFG$cols$cell_type]]
  )

cell_summary <- md %>%
  group_by(donor_id, infection, hpi, condition, cell_type) %>%
  summarise(
    n_cells = dplyr::n(),
    EA_release_score = mean(EA_release_score, na.rm = TRUE),
    EA_recycle_score = mean(EA_recycle_score, na.rm = TRUE),
    EA_proxy = mean(EA_proxy, na.rm = TRUE),
    FLT1_expr = mean(FLT1_expr, na.rm = TRUE),
    PGF_expr = mean(PGF_expr, na.rm = TRUE),
    FLT1_minus_PGF = mean(FLT1_minus_PGF, na.rm = TRUE),
    Angio_stress_score = mean(Angio_stress_score, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(cell_summary, file.path(tbl_dir, "ea_flt1_celltype_condition_summary.csv"), row.names = FALSE)
qs::qsave(cell_summary, file.path(obj_dir, "ea_flt1_celltype_condition_summary.qs"), preset = "high")

fit_mixed_by_celltype <- function(df, response) {
  out <- list()
  for (ct in unique(df$cell_type)) {
    d <- df %>% filter(cell_type == ct)
    if (nrow(d) < 8 || length(unique(d$donor_id)) < 3) next
    d$w <- log1p(d$n_cells)
    form <- stats::as.formula(paste0(response, " ~ infection * hpi + (1|donor_id)"))
    fit <- tryCatch(lmer(form, data = d, weights = w), error = function(e) NULL)
    if (is.null(fit)) next
    an <- suppressMessages(anova(fit))
    an_df <- data.frame(term = rownames(an), an, row.names = NULL)
    an_df$cell_type <- ct
    an_df$response <- response
    out[[ct]] <- an_df
  }
  bind_rows(out)
}

anova_tbl <- bind_rows(
  fit_mixed_by_celltype(cell_summary, "EA_proxy"),
  fit_mixed_by_celltype(cell_summary, "FLT1_minus_PGF"),
  fit_mixed_by_celltype(cell_summary, "Angio_stress_score")
)
write.csv(anova_tbl, file.path(tbl_dir, "ea_flt1_mixed_model_anova.csv"), row.names = FALSE)

heat_df <- cell_summary %>%
  mutate(cond = paste(infection, hpi, sep = "_")) %>%
  select(cell_type, cond, EA_proxy, FLT1_minus_PGF, Angio_stress_score) %>%
  pivot_longer(cols = c(EA_proxy, FLT1_minus_PGF, Angio_stress_score), names_to = "metric", values_to = "value") %>%
  group_by(cell_type, cond, metric) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

p_heat <- ggplot(heat_df, aes(x = cond, y = cell_type, fill = value)) +
  geom_tile() +
  facet_wrap(~ metric, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell-type by condition summary of EA and FLT1-related host proxies", x = "Condition", y = "Cell type")

save_plot(file.path(fig_dir, "Fig27_Heatmap_EA_FLT1_by_celltype_condition.pdf"), p_heat, w = 16, h = 8)

spaghetti_ct <- intersect(c("VCT_fusing", "VCT", "EVT_1", "EVT_2", "F", "PV", "Endo_f", "HBC"), unique(cell_summary$cell_type))
spag_df <- cell_summary %>% filter(cell_type %in% spaghetti_ct)

plot_spaghetti <- function(df, yvar, title_txt, outfile) {
  p <- ggplot(df, aes(x = condition, y = .data[[yvar]], group = donor_id, color = donor_id)) +
    geom_line(alpha = 0.6) +
    geom_point(size = 1.6) +
    facet_wrap(~ cell_type, scales = "free_y") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title_txt, x = "Condition", y = yvar)
  save_plot(outfile, p, w = 16, h = 10)
}

plot_spaghetti(spag_df, "EA_proxy", "Matched-donor trajectories of EA proxy", file.path(fig_dir, "Fig28_Spaghetti_EA_proxy.pdf"))
plot_spaghetti(spag_df, "FLT1_minus_PGF", "Matched-donor trajectories of FLT1 minus PGF proxy", file.path(fig_dir, "Fig29_Spaghetti_FLT1minusPGF.pdf"))
plot_spaghetti(spag_df, "Angio_stress_score", "Matched-donor trajectories of angiogenic stress score", file.path(fig_dir, "Fig30_Spaghetti_AngioStress.pdf"))

qs::qsave(seu, file.path(obj_dir, "seu_with_ea_flt1_proxies.qs"), preset = "high")

write_legend(
  "Fig21-30", "EA-availability and FLT1-axis proxy analyses",
  hypothesis = "Host-side ethanolamine liberation programs and angiogenic stress programs vary by cell type, infection, and time, and these shifts may reveal placental states permissive for an Fn-like metabolic mechanism.",
  methods = "We computed sparse host proxy scores for EA liberation, EA sink, and FLT1-related angiogenic stress, summarized scores by donor, condition, and cell type, and fit mixed-effects models with donor as a random effect.",
  readout = "Preferred signals are consistent cell-type-restricted shifts that remain visible in exact genes, score summaries, and matched-donor trajectories. Total FLT1 is measured here; the placental-specific sFLT1 e15a isoform is not directly resolved by standard scRNA-seq.",
  interpretation_template = "- Strong EA proxy elevation in trophoblast or fibroblast compartments supports increased substrate availability.\n- Increased FLT1-minus-PGF or angiogenic stress in trophoblast/endothelium supports a PE-relevant vascular axis.\n- Signals confined to a few donors are hypothesis-generating, not definitive.",
  outfile = file.path(CFG$dirs$legends, "Fig21_30_legend.md")
)

log_msg("10_ea_sflt1_proxy_analysis done.", log_file = log_file)
