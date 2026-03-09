source("R/config.R")
source("R/helpers_io.R")
source("R/helpers_plot.R")
source("R/helpers_scores.R")

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("06_nk_cytotoxic_module starting.", log_file = log_file)

seu <- qs::qread(file.path(CFG$dirs$objects, "seu_with_scores.qs"))
DefaultAssay(seu) <- "RNA"

# IMPORTANT NOTE:
# `seu$GENE` accesses metadata, not expression. Use FetchData() for expression.

markers <- c("PTPRC","NKG7","KLRD1","GNLY","PRF1","GZMB","FCGR3A",
             "TRAC","CD3D","CD3E","IL7R",
             "LST1","TYROBP","LYZ",
             "KRT8","KRT18")
markers <- present_genes(seu, markers)

p13 <- FeaturePlot(seu, features = markers, cols = c("lightgrey","blue"), order = TRUE, ncol = 4, raster = TRUE)
save_plot(file.path(CFG$dirs$figures, "Fig13_FeaturePlot_cytotoxic_panel.pdf"), p13, w = 16, h = 10)

write_legend(
  "Fig13", "Marker panel for NK/T/myeloid/trophoblast identity",
  hypothesis = "If true NK/T cells are present, they should show coherent immune identity (PTPRC+) and lineage markers (TRAC/CD3D for T; NKG7/KLRD1/GNLY/PRF1 for NK) without trophoblast contamination (KRT8/18).",
  methods = "Overlay immune and cytotoxic markers on the global UMAP; check co-localization and exclusivity.",
  readout = "A real lymphocyte cluster would be PTPRC+ and TRAC/CD3D+ (T) or PTPRC+ and NKG7/KLRD1/GNLY+ (NK), and largely KRT-.",
  interpretation_template = "- If markers are sparse and not co-expressed, likely ambient RNA or doublets.\n- If a coherent cluster exists, subset and re-cluster that compartment for deeper profiling."
)

expr <- FetchData(seu, vars = markers)

seu$cytotoxic_like <- with(expr, (NKG7 > 0 | GNLY > 0 | PRF1 > 0 | GZMB > 0))
seu$immune_cytotoxic_like <- with(expr, PTPRC > 0 & (NKG7 > 0 | GNLY > 0 | PRF1 > 0 | GZMB > 0))

seu$NK_like_strict <- with(expr,
                           PTPRC > 0 &
                             NKG7 > 0 & KLRD1 > 0 & (GNLY > 0 | PRF1 > 0) &
                             TRAC == 0 & CD3D == 0 & CD3E == 0 & IL7R == 0)

# Abundance by condition

tab <- seu@meta.data %>%
  count(infection = .data[[CFG$cols$infection]], hpi = .data[[CFG$cols$hpi]], NK_like_strict) %>%
  group_by(infection, hpi) %>%
  mutate(frac = n/sum(n)) %>% ungroup()

write.csv(tab, file.path(CFG$dirs$tables, "NK_like_strict_by_condition.csv"), row.names = FALSE)

p14 <- DimPlot(seu, group.by = "NK_like_strict") + ggtitle("NK_like_strict flagged cells")
save_plot(file.path(CFG$dirs$figures, "Fig14_UMAP_NK_like_strict.pdf"), p14, w = 12, h = 8)

write_legend(
  "Fig14", "UMAP: NK_like_strict flagged cells",
  hypothesis = "If infection recruits/expands NK-like cells or induces NK-like cytotoxic states, flagged cells should increase by infection/time and form a coherent cluster.",
  methods = "Compute strict NK-like boolean gate from expression fetched via FetchData; plot flagged cells on UMAP; tabulate by condition.",
  readout = "A coherent cluster + condition-associated increase supports a real signal; scattered rare cells suggests contamination/doublets.",
  interpretation_template = "- If extremely rare (<0.1%) and scattered, treat as likely contamination.\n- If enriched in a condition, re-cluster immune cells and validate with additional markers (FCGR3A, TRBC1/2, TRAC, MS4A7, etc.)."
)

p15 <- VlnPlot(seu, features = c("nCount_RNA","nFeature_RNA","percent.mt"), group.by = "NK_like_strict", pt.size = 0) +
  ggtitle("QC metrics for NK_like_strict vs others")
save_plot(file.path(CFG$dirs$figures, "Fig15_QC_NK_like_strict.pdf"), p15, w = 12, h = 5)

write_legend(
  "Fig15", "QC: NK_like_strict vs others",
  hypothesis = "If flagged NK-like cells are doublets, they often show elevated nCount_RNA / nFeature_RNA.",
  methods = "Compare QC distributions for NK_like_strict vs non-flagged cells.",
  readout = "Higher counts/features in flagged group suggests doublets; similar QC supports real cells.",
  interpretation_template = "- If doublet-like, remove or ignore for biological inference.\n- If QC normal and markers coherent, proceed to subset/recluster for activation programs."
)

qs::qsave(seu, file.path(CFG$dirs$objects, "seu_with_nk_flags.qs"))
log_msg("06_nk_cytotoxic_module done.", log_file = log_file)
