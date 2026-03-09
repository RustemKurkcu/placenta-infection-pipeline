source("R/config.R")
source("R/helpers_io.R")
source("R/helpers_plot.R")

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("03_core_umaps_and_composition starting.", log_file = log_file)

seu <- qs::qread(file.path(CFG$dirs$objects, "seu_qc_checked.qs"))
DefaultAssay(seu) <- "RNA"

if (!"umap" %in% names(seu@reductions) && !"X_umap" %in% names(seu@reductions)) {
  log_msg("No UMAP found; running standard normalization + PCA + UMAP.", log_file = log_file)
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
  seu <- RunPCA(seu, features = VariableFeatures(seu), verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
}

p3 <- DimPlot(seu, group.by = CFG$cols$cell_type, label = TRUE, repel = TRUE) +
  ggtitle("UMAP colored by cell type")
save_plot(file.path(CFG$dirs$figures, "Fig03_UMAP_cell_type.pdf"), p3, w = 12, h = 8)

write_legend(
  "Fig03", "UMAP by cell type",
  hypothesis = "Placental compartments (trophoblast, macrophage, stromal, endothelial) segregate into distinct transcriptomic states.",
  methods = "Plot existing UMAP embedding and color by the provided cell_type annotation.",
  readout = "Well-separated, interpretable clusters suggest good annotation and biology; mixing may indicate batch effects or ambiguous states.",
  interpretation_template = "- If cell types are well separated, proceed with within-cell-type comparisons.\n- If not, consider re-clustering or checking marker genes."
)

p4 <- dimplot_split(seu, group_by = CFG$cols$infection, title = "UMAP by infection")
save_plot(file.path(CFG$dirs$figures, "Fig04_UMAP_by_infection.pdf"), p4, w = 12, h = 8)

write_legend(
  "Fig04", "UMAP by infection",
  hypothesis = "Infection may induce global state shifts or introduce condition-specific subclusters.",
  methods = "DimPlot colored by infection label.",
  readout = "Condition-enriched regions may reflect infection-driven transcriptional states or sampling differences.",
  interpretation_template = "- If certain infections occupy distinct regions, follow up within each cell type for DE and pathway signals.\n- If infections overlap, effects may be subtle or cell-type specific."
)

p5 <- dimplot_split(seu, group_by = "hpi_num", title = "UMAP by hours post infection (numeric)")
save_plot(file.path(CFG$dirs$figures, "Fig05_UMAP_by_time.pdf"), p5, w = 12, h = 8)

write_legend(
  "Fig05", "UMAP by time (hpi)",
  hypothesis = "Host responses evolve from 24h to 48h; time may separate early vs late activation states.",
  methods = "Convert hpi (e.g., '24h','48h') to numeric and color UMAP by hpi_num.",
  readout = "Gradients or enriched regions at 48h suggest maturation/activation of responses.",
  interpretation_template = "- If time separates strongly, model time as a key factor and test interaction with infection.\n- If weak, focus on within-cell-type gene programs."
)

# Composition by infection

df_comp <- seu@meta.data %>%
  count(.data[[CFG$cols$infection]], .data[[CFG$cols$cell_type]]) %>%
  group_by(.data[[CFG$cols$infection]]) %>%
  mutate(frac = n / sum(n)) %>% ungroup()

p6 <- ggplot(df_comp, aes(x = .data[[CFG$cols$infection]], y = frac, fill = .data[[CFG$cols$cell_type]])) +
  geom_col(width = 0.8) +
  ylab("Fraction of cells") + xlab("Infection") +
  ggtitle("Cell-type composition by infection") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(file.path(CFG$dirs$figures, "Fig06_Composition_by_infection.pdf"), p6, w = 12, h = 8)

write_legend(
  "Fig06", "Cell-type composition by infection",
  hypothesis = "Some infections may alter survival/representation of specific compartments (immune enrichment, trophoblast depletion) or reflect selective infection of niches.",
  methods = "Count cells per infection × cell_type and convert to fractions per infection.",
  readout = "Large shifts for a given infection suggest differential survival, dissociation bias, or true compositional remodeling.",
  interpretation_template = "- Validate with per-sample (donor) composition to rule out sampling artifacts.\n- If stable, prioritize state changes (gene programs) rather than abundance changes."
)

ensure_dir(CFG$dirs$tables)
write.csv(df_comp, file.path(CFG$dirs$tables, "composition_by_infection.csv"), row.names = FALSE)

qs::qsave(seu, file.path(CFG$dirs$objects, "seu_core_figs.qs"))
log_msg("03_core_umaps_and_composition done.", log_file = log_file)
