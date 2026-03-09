source("R/config.R")
source("R/helpers_io.R")
source("R/helpers_plot.R")
source("R/helpers_scores.R")

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("04_gene_sets_scores_plots starting.", log_file = log_file)

seu <- qs::qread(file.path(CFG$dirs$objects, "seu_core_figs.qs"))
DefaultAssay(seu) <- "RNA"

# Load gene sets

glyco_genes    <- readLines("gene_sets/glyco_genes.txt")
adhesion_genes <- readLines("gene_sets/adhesion_genes.txt")
innate_genes   <- readLines("gene_sets/innate_genes.txt")

# keep only present genes

glyco_genes    <- present_genes(seu, glyco_genes)
adhesion_genes <- present_genes(seu, adhesion_genes)
innate_genes   <- present_genes(seu, innate_genes)

# Dotplots
p7 <- DotPlot(seu, features = glyco_genes, group.by = CFG$cols$cell_type) + RotatedAxis() +
  ggtitle("Glycocalyx / O-glyco proxy genes by cell type")
save_plot(file.path(CFG$dirs$figures, "Fig07_Dotplot_glyco_by_celltype.pdf"), p7, w = 14, h = 6)

write_legend(
  "Fig07", "DotPlot: glycocalyx / O-glycosylation proxies",
  hypothesis = "Barrier/entry-relevant glycocalyx programs are enriched in specific trophoblast compartments (candidate susceptible niches).",
  methods = "DotPlot summarizing average expression (color) and percent expressing (dot size) of curated glyco/glycocalyx genes per cell type.",
  readout = "Cell types with high mucin/proteoglycan/glycosyltransferase signal.",
  interpretation_template = "- Enrichment in trophoblast supports a susceptibility proxy.\n- If mostly immune/stromal, refine the gene set or score within trophoblast only."
)

p8 <- DotPlot(seu, features = adhesion_genes, group.by = CFG$cols$cell_type) + RotatedAxis() +
  ggtitle("Adhesion / epithelial interaction proxy genes by cell type")
save_plot(file.path(CFG$dirs$figures, "Fig08_Dotplot_adhesion_by_celltype.pdf"), p8, w = 11, h = 5)

write_legend(
  "Fig08", "DotPlot: adhesion proxies",
  hypothesis = "Epithelial adhesion programs (CDH1/EPCAM/integrins) mark potential pathogen-binding/entry-permissive trophoblast states.",
  methods = "DotPlot of curated adhesion genes by cell type.",
  readout = "Cell types with highest CDH1/EPCAM/integrin signal.",
  interpretation_template = "- High adhesion in trophoblast supports entry-susceptibility framing.\n- If uniform, use pathogen-specific receptor sets."
)

p9 <- DotPlot(seu, features = innate_genes, group.by = CFG$cols$cell_type) + RotatedAxis() +
  ggtitle("Innate / inflammatory proxy genes by cell type")
save_plot(file.path(CFG$dirs$figures, "Fig09_Dotplot_innate_by_celltype.pdf"), p9, w = 11, h = 5)

write_legend(
  "Fig09", "DotPlot: innate/inflammatory proxies",
  hypothesis = "Infection severity manifests as elevated innate programs in immune compartments and/or stressed trophoblast.",
  methods = "DotPlot of inflammatory and interferon-response proxy genes by cell type.",
  readout = "Cell types with strongest innate baseline/inducible signal.",
  interpretation_template = "- Strong signal in macrophage-like cells supports an immune severity readout.\n- If trophoblast shows strong induction, severity may be epithelial-intrinsic."
)

# Signature scores (memory safe)
seu <- add_signature_hvg_pool(seu, glyco_genes, "GlycoScore")
seu <- add_signature_hvg_pool(seu, adhesion_genes, "AdhesionScore")
seu <- add_signature_hvg_pool(seu, innate_genes, "InnateScore")

# Optional within-cell-type z scores
seu <- zscore_within(seu, "GlycoScore",    CFG$cols$cell_type)
seu <- zscore_within(seu, "AdhesionScore", CFG$cols$cell_type)
seu <- zscore_within(seu, "InnateScore",   CFG$cols$cell_type)

# White->red FeaturePlot for clarity
p10 <- featureplot_white_red(seu, "AdhesionScore", title = "UMAP: AdhesionScore (white→red)")
save_plot(file.path(CFG$dirs$figures, "Fig10_UMAP_AdhesionScore_white_red.pdf"), p10, w = 12, h = 8)

write_legend(
  "Fig10", "UMAP overlay: AdhesionScore",
  hypothesis = "Entry/attachment susceptibility varies across cell states; adhesion-high niches may be more permissive.",
  methods = "Compute an adhesion signature per cell and overlay on UMAP; use quantile cutoffs and ordering for clarity.",
  readout = "Localized high-score regions (red) indicate putative entry-permissive states.",
  interpretation_template = "- If adhesion-high aligns with specific trophoblast subtypes, prioritize them for infection/DE analysis.\n- If widespread, adhesion is not discriminatory; refine gene set."
)

# Split by infection and time
p11 <- featureplot_white_red(seu, "AdhesionScore", split_by = "condition",
                            title = "AdhesionScore by condition (infection × time)")
save_plot(file.path(CFG$dirs$figures, "Fig11_AdhesionScore_split_condition.pdf"), p11, w = 16, h = 10)

write_legend(
  "Fig11", "AdhesionScore split by condition",
  hypothesis = "Infection/time may shift adhesion programs (either up-regulation or selective survival of adhesion-high cells).",
  methods = "Split FeaturePlot panels by condition (infection × hpi).",
  readout = "Panel-to-panel shifts in intensity or spatial localization.",
  interpretation_template = "- If adhesion increases in specific infections over time, consider it a response program.\n- If only certain infections show a shift, investigate pathogen-specific entry mechanisms."
)

qs::qsave(seu, file.path(CFG$dirs$objects, "seu_with_scores.qs"))
log_msg("04_gene_sets_scores_plots done.", log_file = log_file)
