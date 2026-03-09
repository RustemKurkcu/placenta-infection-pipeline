# ==============================================================================
# MASTER SCRIPT: PLACENTA INFECTION PIPELINE (Explant -> Organoid Prep)
# ==============================================================================
# MERGES: User's 'placenta_fn_pipeline.R' + Fixes for tibble/edgeR + Next Steps
# OUTPUTS: Figures, Tables, Legends (Markdown), and Logs
# ==============================================================================

# 1. SETUP & LIBRARIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(edgeR)
  library(qs)
  library(ggrepel) # For Volcano plots
})

# --- Configuration ---
CFG <- list(
  qs_seurat = "data/02_processed/seurat_object.qs",
  fig_dir   = "results/figures",
  tab_dir   = "results/tables",
  leg_dir   = "results/legends",
  log_dir   = "results/logs"
)

# --- Create Dirs ---
for(d in c(CFG$fig_dir, CFG$tab_dir, CFG$leg_dir, CFG$log_dir)) {
  dir.create(d, recursive=TRUE, showWarnings=FALSE)
}

# --- Helper Functions ---
log_msg <- function(...) {
  cat(format(Sys.time(), "%F %T"), "|", paste(..., collapse=" "), "\n",
      file = file.path(CFG$log_dir, "run.log"), append = TRUE)
  message(paste(..., collapse=" "))
}

save_plot <- function(filename, plot_obj, w=8, h=6) {
  ggsave(file.path(CFG$fig_dir, filename), plot_obj, width=w, height=h)
  log_msg("Saved Plot:", filename)
}

write_legend <- function(fig_id, title, hypothesis, methods, readout) {
  txt <- paste0(
    "# ", fig_id, ": ", title, "\n\n",
    "## Hypothesis\n", hypothesis, "\n\n",
    "## Methods\n", methods, "\n\n",
    "## Readout\n", readout, "\n"
  )
  writeLines(txt, file.path(CFG$leg_dir, paste0(fig_id, ".md")))
}

# 2. LOAD DATA
# ------------------------------------------------------------------------------
log_msg("--- STARTING PIPELINE ---")
if (file.exists(CFG$qs_seurat)) {
  seu <- qs::qread(CFG$qs_seurat)
  DefaultAssay(seu) <- "RNA"
  log_msg("Loaded Seurat Object with", ncol(seu), "cells.")
} else {
  stop("Seurat object not found at ", CFG$qs_seurat)
}

# 3. CORE FIGURES (QC & UMAP)
# ------------------------------------------------------------------------------
# Fig 01: QC
p_qc <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
save_plot("Fig01_QC.pdf", p_qc, w=10, h=5)
write_legend("Fig01", "Quality Control", "Data is clean enough for analysis.", "Violin plots of metrics.", "No massive dropout.")

# Fig 02: UMAP Cell Type
p_umap <- DimPlot(seu, group.by = "cell_type", label = TRUE) + ggtitle("Cell Types")
save_plot("Fig02_UMAP_CellType.pdf", p_umap, w=10, h=8)
write_legend("Fig02", "Cell Type UMAP", "Distinct lineages exist.", "Standard Seurat UMAP.", "Clear clusters.")

# 4. MODULE SCORING (Susceptibility vs Severity)
# ------------------------------------------------------------------------------
log_msg("Calculating Module Scores...")
genes_susc <- intersect(c("CDH1", "EPCAM", "ITGA6", "MUC1", "SDC1", "GALNT1"), rownames(seu))
genes_sev  <- intersect(c("IL1B", "TNF", "CXCL8", "CCL3", "NFKBIA", "ISG15"), rownames(seu))

seu <- AddModuleScore(seu, features = list(genes_susc), name = "SusceptibilityScore")
seu <- AddModuleScore(seu, features = list(genes_sev), name = "SeverityScore")
seu$Susc_Score <- seu$SusceptibilityScore1
seu$Sev_Score  <- seu$SeverityScore1

# Fig 09: The Correlation Plot
# We aggregate by Sample (Donor + Stage) to see the trend
df_susc <- seu@meta.data %>% filter(cell_type %in% c("VCT", "EVT_1", "SCT")) %>%
  group_by(sample_id = paste(donor_id, stage, sep="_")) %>% summarise(Mean_Susc = mean(Susc_Score, na.rm=TRUE))

df_sev <- seu@meta.data %>% filter(cell_type %in% c("HBC", "PAMM1")) %>%
  group_by(sample_id = paste(donor_id, stage, sep="_")) %>% summarise(Mean_Sev = mean(Sev_Score, na.rm=TRUE))

plot_data <- inner_join(df_susc, df_sev, by="sample_id") %>%
  mutate(Infection = ifelse(grepl("UI", sample_id), "UI", "Infected"))

p_corr <- ggplot(plot_data, aes(x=Mean_Susc, y=Mean_Sev, color=Infection)) +
  geom_point(size=3) + geom_smooth(method="lm", se=FALSE) + theme_bw() +
  labs(title="Susceptibility vs Severity", subtitle="Hypothesis Test")

save_plot("Fig09_Correlation.pdf", p_corr)
write_legend("Fig09", "Susceptibility vs Severity", 
             "Higher adhesion markers in Trophs lead to higher inflammation in HBCs.",
             "Scatterplot of module scores aggregated by sample.",
             "Positive correlation would support hypothesis. Flat/Negative refutes it.")

# 5. PSEUDOBULK DE (The Fixed Logic)
# ------------------------------------------------------------------------------
log_msg("Starting Pseudobulk DE...")

target_cell_types <- c("VCT", "EVT_1", "HBC", "PAMM1")
contrasts <- c("Lm_24h - UI_24h", "Lm_48h - UI_48h") # Add Pf/Tg if needed

run_de <- function(seu_obj, cell_type) {
  log_msg("Processing DE for:", cell_type)
  sub <- subset(seu_obj, subset = cell_type == cell_type)
  if (ncol(sub) < 50) return(NULL)
  
  # 1. Aggregate
  sub$sample_id <- paste(sub$donor_id, sub$stage, sep = "_")
  cts <- AggregateExpression(sub, group.by = "sample_id", assays = "RNA", slot = "counts")$RNA
  
  # 2. Meta (FIXED TIBBLE ERROR HERE)
  meta <- sub@meta.data %>% distinct(sample_id, stage, donor_id) %>%
    tibble::remove_rownames() %>% tibble::column_to_rownames("sample_id")
  meta <- meta[colnames(cts), ]
  
  # 3. edgeR
  y <- DGEList(counts=cts, group=meta$stage)
  keep <- filterByExpr(y, group=meta$stage)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  
  design <- model.matrix(~ 0 + stage + donor_id, data=meta)
  colnames(design) <- gsub("stage", "", colnames(design))
  colnames(design) <- gsub("donor_id", "donor", colnames(design))
  
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  
  # 4. Contrasts
  my_contrasts <- makeContrasts(contrasts=contrasts, levels=design)
  
  for(k in colnames(my_contrasts)) {
    try({
      qlf <- glmQLFTest(fit, contrast=my_contrasts[,k])
      res <- topTags(qlf, n=Inf)$table
      
      # Save Table
      fname <- paste0("DE_", cell_type, "_", k, ".csv")
      write.csv(res, file.path(CFG$tab_dir, fname))
      log_msg("Saved DE Table:", fname)
      
      # 6. NEXT STEPS: VOLCANO PLOTS (New Feature)
      # -------------------------------------------------
      res$gene <- rownames(res)
      res$Sig <- ifelse(res$FDR < 0.05 & abs(res$logFC) > 1, "Significant", "NS")
      
      p_vol <- ggplot(res, aes(x=logFC, y=-log10(PValue), color=Sig)) +
        geom_point(alpha=0.5) + theme_minimal() +
        scale_color_manual(values=c("grey", "red")) +
        geom_text_repel(data=head(res, 10), aes(label=gene), color="black") +
        ggtitle(paste("Volcano:", cell_type, k))
      
      save_plot(paste0("Fig_Volcano_", cell_type, "_", k, ".pdf"), p_vol)
      
    }, silent=TRUE)
  }
}

for(ct in target_cell_types) run_de(seu, ct)

# 7. NK CELL CHECK (From your script)
# ------------------------------------------------------------------------------
log_msg("Checking for NK cells...")
genes_nk <- c("PTPRC","NKG7","KLRD1","GNLY","PRF1")
if(all(genes_nk %in% rownames(seu))) {
  seu$NK_Score <- colMeans(GetAssayData(seu)[genes_nk,])
  p_nk <- FeaturePlot(seu, features="NK_Score") + ggtitle("NK Cell Markers")
  save_plot("Fig_NK_Check.pdf", p_nk)
}

log_msg("--- PIPELINE COMPLETE ---")