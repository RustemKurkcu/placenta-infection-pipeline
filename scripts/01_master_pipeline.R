# ==============================================================================
# 01_master_pipeline.R
# Placenta infection scRNA-seq reanalysis (Hoo et al., Cell Systems 2024) -> Fn organoid test-case
#
# What this script does (high level):
#   1) Loads a processed Seurat object (qs) with metadata already attached
#   2) Generates QC + UMAP + composition figures with saved legends
#   3) Computes "proxy" module scores in a memory-safe (sparse) way
#   4) Tests the "Susceptibility -> Severity" hypothesis at the SAMPLE level
#   5) Runs per-cell-type pseudobulk DE with edgeR using *matched* uninfected controls
#
# Key design choice:
#   Prefer stage_perInfection (e.g., UI_Tg_24h) over stage (UI_24h) so that each pathogen
#   is compared to its matched control prepared in parallel.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(Matrix)
  library(edgeR)
  library(qs)
  library(ggrepel)
})

# -------------------------------
# CONFIG
# -------------------------------
CFG <- list(
  qs_seurat = "data/02_processed/seurat_object.qs",
  out_dir   = "results",
  fig_dir   = "results/figures",
  tab_dir   = "results/tables",
  leg_dir   = "results/legends",
  log_dir   = "results/logs"
)

dir.create(CFG$out_dir, recursive=TRUE, showWarnings=FALSE)
for(d in c(CFG$fig_dir, CFG$tab_dir, CFG$leg_dir, CFG$log_dir)) dir.create(d, recursive=TRUE, showWarnings=FALSE)

log_msg <- function(...) {
  msg <- paste(format(Sys.time(), "%F %T"), "|", paste(..., collapse=" "))
  cat(msg, "\n", file=file.path(CFG$log_dir, "run.log"), append=TRUE)
  message(msg)
}

save_plot <- function(filename, plot_obj, w=8, h=6) {
  ggsave(filename = file.path(CFG$fig_dir, filename), plot = plot_obj, width=w, height=h, useDingbats=FALSE)
  log_msg("Saved plot:", filename)
}

write_legend <- function(fig_id, title, hypothesis, methods, readout, interpretation=NULL) {
  txt <- paste0(
    "# ", fig_id, ": ", title, "\n\n",
    "## Hypothesis\n", hypothesis, "\n\n",
    "## Methods\n", methods, "\n\n",
    "## Readout\n", readout, "\n\n",
    if (!is.null(interpretation)) paste0("## How to interpret\n", interpretation, "\n") else ""
  )
  writeLines(txt, con=file.path(CFG$leg_dir, paste0(fig_id, ".md")))
  log_msg("Saved legend:", paste0(fig_id, ".md"))
}

# -------------------------------
# LOAD DATA
# -------------------------------
log_msg("=== START ===")
stopifnot(file.exists(CFG$qs_seurat))
seu <- qs::qread(CFG$qs_seurat)
DefaultAssay(seu) <- "RNA"

log_msg("Loaded Seurat object:", ncol(seu), "cells,", nrow(seu), "genes")

# Choose grouping variable for matched controls
group_var <- if ("stage_perInfection" %in% colnames(seu@meta.data)) "stage_perInfection" else "stage"
log_msg("Using group variable:", group_var)

# Make sure basic fields exist
stopifnot(all(c("cell_type","infection","hpi","donor_id") %in% colnames(seu@meta.data)))

# -------------------------------
# FIG 01: QC
# -------------------------------
p_qc <- VlnPlot(seu, features=c("nCount_RNA","nFeature_RNA","percent.mt"), ncol=3, pt.size=0)
save_plot("Fig01_QC.pdf", p_qc, w=12, h=4)
write_legend(
  "Fig01", "QC metrics",
  hypothesis="Basic QC metrics are within expected ranges; no sample is dominated by stressed/dying cells or extreme doublet-like libraries.",
  methods="Violin plots of nCount_RNA, nFeature_RNA, and percent.mt across all cells.",
  readout="Broadly similar distributions; absence of extreme high percent.mt tails; no pervasive ultra-high nCount outliers.",
  interpretation="If any group shows extreme percent.mt or nCount/nFeature tails, apply cautious filtering and re-run."
)

# -------------------------------
# FIG 02: UMAP by cell type
# -------------------------------
p_umap_ct <- DimPlot(seu, group.by="cell_type", label=TRUE, repel=TRUE) + ggtitle("UMAP: cell types")
save_plot("Fig02_UMAP_celltype.pdf", p_umap_ct, w=11, h=8)
write_legend(
  "Fig02", "UMAP by cell type",
  hypothesis="Major placental lineages (trophoblast, macrophage) form coherent clusters, enabling within-lineage comparisons.",
  methods="DimPlot on existing embedding, colored by cell_type.",
  readout="Distinct clusters for VCT/EVT and HBC/PAMM with minimal mixing.",
  interpretation="Good separation supports downstream DE within cell types (state changes rather than re-clustering artifacts)."
)

# -------------------------------
# FIG 03: UMAP by infection/time (group_var)
# -------------------------------
p_umap_stage <- DimPlot(seu, group.by=group_var, shuffle=TRUE) + ggtitle(paste0("UMAP: ", group_var))
save_plot("Fig03_UMAP_stage.pdf", p_umap_stage, w=12, h=8)
write_legend(
  "Fig03", paste0("UMAP by ", group_var),
  hypothesis="Infection does not globally rewrite cell identity; changes are expected as gene-program shifts within lineages.",
  methods="DimPlot colored by experimental condition.",
  readout="Conditions are intermingled within clusters (no condition-specific 'new cluster' dominating the map).",
  interpretation="If one condition forms a distinct island, check QC, doublets, or strong cell-state transitions."
)

# -------------------------------
# FIG 04: Composition (per sample_id, then per condition)
# -------------------------------
# Build a robust sample_id if missing
if (!"sample_id" %in% colnames(seu@meta.data)) {
  seu$sample_id <- paste(seu$donor_id, seu[[group_var, drop=TRUE]], sep="|")
}

df_comp <- seu@meta.data %>%
  count(sample_id, donor_id, !!sym(group_var), infection, hpi, cell_type) %>%
  group_by(sample_id) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

p_comp <- ggplot(df_comp, aes(x=cell_type, y=frac, fill=infection)) +
  geom_col(position="stack") +
  facet_wrap(~ !!sym(group_var), ncol=4) +
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(title="Cell-type composition by condition", y="Fraction of cells", x="Cell type")

save_plot("Fig04_composition_by_condition.pdf", p_comp, w=16, h=9)
write_legend(
  "Fig04", "Cell-type composition by condition",
  hypothesis="Pathogens might alter cell representation (loss/enrichment) or mostly alter state within stable composition.",
  methods="Counts per sample_id × cell_type converted to fractions; faceted by condition.",
  readout="Large shifts in a lineage fraction suggest compositional remodeling or technical bias.",
  interpretation="If composition is stable, prioritize within-cell-type DE and signaling analyses."
)

# -------------------------------
# PROXY GENE SETS (edit as needed)
# -------------------------------
present_genes <- function(seu_obj, genes) intersect(genes, rownames(seu_obj))

glyco_genes <- present_genes(seu, c(
  "GALNT1","GALNT2","GALNT3","GALNT7","GALNT10",
  "C1GALT1","C1GALT1C1","ST3GAL1","ST6GALNAC1","ST6GALNAC2",
  "MUC1","MUC4","SDC1","SDC4","HSPG2","VCAN"
))

adhesion_genes <- present_genes(seu, c("CDH1","EPCAM","ITGA6","ITGB1","ICAM1","VCAM1"))
innate_genes   <- present_genes(seu, c("TLR2","TLR4","NFKBIA","IL1B","CXCL8","TNF","IL6","ISG15","IFIT1","IFIT2","IFIT3","RSAD2","MX1"))

# Dotplots (cell-type context)
p_glyco <- DotPlot(seu, features=glyco_genes, group.by="cell_type") + RotatedAxis() +
  ggtitle("Glycocalyx / O-glyco proxy genes by cell_type")
save_plot("Fig05_dotplot_glyco_by_celltype.pdf", p_glyco, w=16, h=6)

p_adh <- DotPlot(seu, features=adhesion_genes, group.by="cell_type") + RotatedAxis() +
  ggtitle("Adhesion proxy genes by cell_type")
save_plot("Fig06_dotplot_adhesion_by_celltype.pdf", p_adh, w=12, h=5)

p_inn <- DotPlot(seu, features=innate_genes, group.by="cell_type") + RotatedAxis() +
  ggtitle("Innate / IFN proxy genes by cell_type")
save_plot("Fig07_dotplot_innate_by_celltype.pdf", p_inn, w=14, h=6)

write_legend(
  "Fig05", "DotPlot: glyco proxies by cell type",
  hypothesis="Entry-relevant glycocalyx programs are enriched in specific trophoblast compartments (candidate ‘susceptible’ niches).",
  methods="DotPlot summarizes mean expression (color) and percent expressing (dot size) per cell type for curated glycocalyx/glyco genes.",
  readout="Which trophoblast subtypes show the strongest mucin/proteoglycan/glycosyltransferase signal."
)

write_legend(
  "Fig06", "DotPlot: adhesion proxies by cell type",
  hypothesis="Adhesion/epithelial interaction programs mark entry-permissive trophoblast states.",
  methods="DotPlot of curated adhesion genes (CDH1/EPCAM/integrins/ICAM/VCAM) by cell type.",
  readout="Which trophoblast subtypes show strongest CDH1/EPCAM/integrin expression."
)

write_legend(
  "Fig07", "DotPlot: innate/IFN proxies by cell type",
  hypothesis="Severity of response manifests as innate/IFN programs in immune compartments and/or stressed trophoblast.",
  methods="DotPlot of inflammatory and interferon-response proxies by cell type.",
  readout="Which cell types carry strongest innate/IFN signal."
)

# -------------------------------
# MEMORY-SAFE MODULE SCORES (sparse)
# -------------------------------
module_score_sparse <- function(seu_obj, genes, slot="data") {
  genes <- intersect(genes, rownames(seu_obj))
  if (length(genes) < 3) return(rep(NA_real_, ncol(seu_obj)))
  mat <- GetAssayData(seu_obj, assay=DefaultAssay(seu_obj), slot=slot)
  Matrix::colMeans(mat[genes, , drop=FALSE])
}

log_msg("Computing sparse module scores...")
seu$GlycoScore1    <- module_score_sparse(seu, glyco_genes)
seu$AdhesionScore1 <- module_score_sparse(seu, adhesion_genes)
seu$InnateScore1   <- module_score_sparse(seu, innate_genes)

# -------------------------------
# FIG 08: Feature plots for module scores (white->red)
# -------------------------------
p_adh_fp <- FeaturePlot(
  seu, features="AdhesionScore1",
  cols=c("white","red"),
  min.cutoff="q05", max.cutoff="q95",
  order=TRUE, raster=TRUE
) + ggtitle("UMAP: AdhesionScore1 (white→red)")
save_plot("Fig08_FeaturePlot_AdhesionScore.pdf", p_adh_fp, w=10, h=8)

p_inn_fp <- FeaturePlot(
  seu, features="InnateScore1",
  cols=c("white","red"),
  min.cutoff="q05", max.cutoff="q95",
  order=TRUE, raster=TRUE
) + ggtitle("UMAP: InnateScore1 (white→red)")
save_plot("Fig08b_FeaturePlot_InnateScore.pdf", p_inn_fp, w=10, h=8)

write_legend(
  "Fig08", "FeaturePlot: AdhesionScore1",
  hypothesis="Adhesion/entry proxy is spatially enriched in trophoblast regions of the embedding.",
  methods="Sparse module score = mean normalized expression of curated adhesion genes; plotted on UMAP with white→red gradient (clipped to q05–q95).",
  readout="High scores should localize to trophoblast clusters; diffuse high signal suggests contamination or an overly broad gene set."
)

# Optional: split by condition for clarity (many panels)
p_adh_split <- FeaturePlot(
  seu, features="AdhesionScore1",
  split.by=group_var, ncol=4,
  cols=c("white","red"),
  min.cutoff="q05", max.cutoff="q95",
  order=TRUE, raster=TRUE
) + ggtitle(paste0("AdhesionScore1 split by ", group_var))
save_plot("Fig08c_FeaturePlot_AdhesionScore_split.pdf", p_adh_split, w=16, h=10)

# -------------------------------
# FIG 09: Susceptibility vs Severity (sample-level)
# -------------------------------
# Define compartments
troph_types <- c("VCT","VCT_fusing","VCT_p","VCT_CCC","EVT_1","EVT_2","iEVT")
immune_types <- c("HBC","HBC_p","PAMM1")

df_troph <- seu@meta.data %>%
  filter(cell_type %in% troph_types) %>%
  group_by(sample_id, donor_id, infection, hpi, !!sym(group_var)) %>%
  summarise(Susc = mean(AdhesionScore1, na.rm=TRUE), .groups="drop")

df_imm <- seu@meta.data %>%
  filter(cell_type %in% immune_types) %>%
  group_by(sample_id, donor_id, infection, hpi, !!sym(group_var)) %>%
  summarise(Sev = mean(InnateScore1, na.rm=TRUE), .groups="drop")

df_corr <- inner_join(df_troph, df_imm, by=c("sample_id","donor_id","infection","hpi",group_var))

p_corr <- ggplot(df_corr, aes(x=Susc, y=Sev, color=infection)) +
  geom_point(size=3) +
  geom_smooth(method="lm", se=FALSE) +
  theme_bw() +
  labs(title="Sample-level test: trophoblast Susceptibility vs immune Severity",
       subtitle=paste0("Grouping=", group_var),
       x="Mean trophoblast AdhesionScore1",
       y="Mean immune InnateScore1")

save_plot("Fig09_Susc_vs_Sev.pdf", p_corr, w=10, h=7)

write_legend(
  "Fig09", "Sample-level Susceptibility vs Severity",
  hypothesis="Samples with higher trophoblast entry/susceptibility proxies mount stronger innate/IFN severity in macrophage compartments.",
  methods="Aggregate module scores per sample_id within trophoblast vs immune compartments; regress Severity ~ Susceptibility.",
  readout="Positive slope supports the hypothesis; flat/negative suggests severity is not limited by entry proxies (or is pathogen-specific/nonlinear).",
  interpretation="If overall correlation is weak, re-test within each pathogen/timepoint and consider nonlinear/threshold behavior."
)

# -------------------------------
# NK / Cytotoxic check (do NOT hard-gate by gene existence in metadata)
# -------------------------------
log_msg("Cytotoxic/NK-like check (score-based).")
cyto_genes <- present_genes(seu, c("NKG7","KLRD1","GNLY","PRF1","GZMB","FCGR3A","XCL1","XCL2","TRAC","TRBC1","TRBC2"))
seu$CytotoxicScore1 <- module_score_sparse(seu, cyto_genes)

p_cyto <- FeaturePlot(seu, "CytotoxicScore1", cols=c("white","red"), min.cutoff="q05", max.cutoff="q99", order=TRUE, raster=TRUE) +
  ggtitle("Cytotoxic/NK-like score (broad markers)")
save_plot("Fig10_CytotoxicScore.pdf", p_cyto, w=10, h=8)

write_legend(
  "Fig10", "Cytotoxic/NK-like score",
  hypothesis="If NK/T-like cells exist, cytotoxic markers should co-localize with PTPRC+ immune regions and form a small distinct population.",
  methods="Sparse module score across cytotoxic/T/NK markers; plotted on UMAP.",
  readout="A true NK population should show PTPRC+, NKG7/KLRD1/GNLY/PRF1 signal with low LST1/TYROBP (not myeloid) and low KRT8/18 (not trophoblast)."
)

# -------------------------------
# PSEUDOBULK DE (edgeR)
# -------------------------------
make_contrasts_map <- function(levels_vec) {
  # Expected naming:
  #   infected:  Lm_24h, Pf_24h, Tg_24h (etc)
  #   matched UI: UI_Lm_24h, UI_Pf_24h, UI_Tg_24h (preferred)
  # Fallback UI: UI_24h / UI_48h
  pathogens <- c("Lm","Pf","Tg")
  times <- c("24h","48h")

  cm <- list()
  for(p in pathogens) {
    for(t in times) {
      inf <- paste0(p, "_", t)
      ui_match <- paste0("UI_", p, "_", t)
      ui_fallback <- paste0("UI_", t)

      if (inf %in% levels_vec) {
        if (ui_match %in% levels_vec) {
          cm[[paste0(p, "_", t, "_vs_UI_", p, "_", t)]] <- c(inf, ui_match)
        } else if (ui_fallback %in% levels_vec) {
          cm[[paste0(p, "_", t, "_vs_UI_", t)]] <- c(inf, ui_fallback)
        }
      }
    }
  }
  # within-pathogen time
  for(p in pathogens) {
    if (all(c(paste0(p,"_24h"), paste0(p,"_48h")) %in% levels_vec)) {
      cm[[paste0(p, "_48h_vs_", p, "_24h")]] <- c(paste0(p,"_48h"), paste0(p,"_24h"))
    }
  }
  cm
}

run_pseudobulk_edger <- function(seu_obj, celltype_label, group_var, min_cells=100) {
  log_msg("DE: starting cell type:", celltype_label)

  sub_obj <- subset(seu_obj, subset = cell_type == celltype_label)
  if (ncol(sub_obj) < min_cells) {
    log_msg("  DE skip: too few cells:", ncol(sub_obj))
    return(invisible(NULL))
  }

  # sample_id must uniquely represent a biological replicate
  # Use donor_id + condition (group_var)
  sub_obj$pb_id <- paste(sub_obj$donor_id, sub_obj[[group_var, drop=TRUE]], sep="|")

  agg <- AggregateExpression(sub_obj, group.by="pb_id", assays="RNA", slot="counts")$RNA

  meta <- sub_obj@meta.data %>%
    distinct(pb_id, donor_id, condition = !!sym(group_var), infection, hpi)

  # align meta to counts
  meta <- meta[match(colnames(agg), meta$pb_id), ]
  stopifnot(all(meta$pb_id == colnames(agg)))

  meta$condition <- factor(meta$condition)

  contrasts_map <- make_contrasts_map(levels(meta$condition))
  if (length(contrasts_map) == 0) {
    log_msg("  DE skip: no valid contrasts found for conditions:", paste(levels(meta$condition), collapse=", "))
    return(invisible(NULL))
  }

  y <- DGEList(counts=agg)
  keep <- filterByExpr(y, group=meta$condition)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)

  # donor blocking helps if donors appear in multiple conditions; if not, it is harmless
  design <- model.matrix(~ 0 + condition + donor_id, data=meta)
  colnames(design) <- make.names(colnames(design))

  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust=TRUE)

  for (nm in names(contrasts_map)) {
    pair <- contrasts_map[[nm]]
    # makeContrasts needs the *design column names* (make.names applied)
    c1 <- make.names(paste0("condition", pair[1]))
    c0 <- make.names(paste0("condition", pair[2]))

    if (!(c1 %in% colnames(design) && c0 %in% colnames(design))) {
      log_msg("  Contrast skip (missing design cols):", nm)
      next
    }

    contr <- rep(0, ncol(design)); names(contr) <- colnames(design)
    contr[c1] <- 1
    contr[c0] <- -1

    qlf <- glmQLFTest(fit, contrast=contr)
    res <- topTags(qlf, n=Inf)$table
    res$Gene <- rownames(res)

    out_csv <- paste0("DE_", celltype_label, "_", nm, ".csv")
    write.csv(res, file=file.path(CFG$tab_dir, out_csv), row.names=FALSE)
    log_msg("  Saved DE:", out_csv)

    # Volcano
    res$Sig <- ifelse(res$FDR < 0.05 & abs(res$logFC) > 1, "Significant", "NS")
    p_vol <- ggplot(res, aes(x=logFC, y=-log10(PValue), color=Sig)) +
      geom_point(alpha=0.5) +
      theme_minimal() +
      geom_vline(xintercept=c(-1,1), linetype="dashed", alpha=0.4) +
      geom_hline(yintercept=-log10(0.05), linetype="dashed", alpha=0.4) +
      geom_text_repel(
        data = res %>% arrange(FDR) %>% filter(Sig=="Significant") %>% head(12),
        aes(label=Gene), max.overlaps=12, size=3
      ) +
      labs(title=paste0("Volcano: ", celltype_label, " | ", nm))

    save_plot(paste0("Volcano_", celltype_label, "_", nm, ".pdf"), p_vol, w=8, h=6)
  }

  invisible(TRUE)
}

log_msg("Running pseudobulk DE (edgeR QL) ...")
target_cell_types <- intersect(c("VCT","EVT_1","EVT_2","HBC","PAMM1"), unique(seu$cell_type))
for (ct in target_cell_types) {
  run_pseudobulk_edger(seu, ct, group_var=group_var, min_cells=200)
}

log_msg("=== COMPLETE ===")
