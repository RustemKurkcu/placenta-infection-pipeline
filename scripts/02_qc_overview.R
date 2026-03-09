source("R/config.R")
source("R/helpers_io.R")

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("02_qc_overview starting.", log_file = log_file)

seu <- qs::qread(file.path(CFG$dirs$objects, "seu_clean.qs"))
DefaultAssay(seu) <- "RNA"

p1 <- VlnPlot(seu, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
              ncol = 3, pt.size = 0.0) + ggtitle("QC metrics (all cells)")

save_plot(file.path(CFG$dirs$figures, "Fig01_QC_violin_allcells.pdf"), p1, w = 12, h = 4)

write_legend(
  "Fig01", "QC metrics violin plots",
  hypothesis = "We need to verify data quality and identify outliers that could bias downstream analysis (e.g., doublets, dying cells).",
  methods = "Plot nCount_RNA, nFeature_RNA, and percent.mt across all cells (no filtering yet).",
  readout = "Extremely high nCount/nFeature suggests doublets; high percent.mt suggests stressed/dying cells.",
  interpretation_template = "- If strong tails exist, define filtering thresholds and re-run.\n- If QC is reasonable, proceed without aggressive filtering (avoid removing true biology)."
)

p2 <- VlnPlot(seu, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
              group.by = "condition", ncol = 3, pt.size = 0.0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("QC metrics by condition")

save_plot(file.path(CFG$dirs$figures, "Fig02_QC_violin_by_condition.pdf"), p2, w = 14, h = 6)

write_legend(
  "Fig02", "QC by condition",
  hypothesis = "Some infections/time points may alter RNA content or stress signatures; we want to rule out technical artifacts.",
  methods = "Violin plots of QC metrics grouped by condition (infection × time).",
  readout = "Systematic shifts in QC can confound differential expression and signature scoring.",
  interpretation_template = "- If one condition has much higher counts/mito, consider per-condition normalization or sensitivity analyses.\n- If stable, comparisons are more trustworthy."
)

qs::qsave(seu, file.path(CFG$dirs$objects, "seu_qc_checked.qs"))
log_msg("02_qc_overview done.", log_file = log_file)
