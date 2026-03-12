source("R/config.R")
source("R/helpers_io.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(qs)
  library(ggplot2)
})

# ============================================================
# 15_targeted_cellchat_axes.R
#
# HYPOTHESIS
# ----------
# The Fn hypothesis is strengthened if infection-relevant cell-cell
# communication axes can be quantified in the infection object rather than
# inferred only from literature.
#
# TARGET AXES
# -----------
# - FIB2/F-like stroma -> trophoblast: AREG/HGF -> EGFR/MET
# - FIB2/F-like stroma -> endothelium: VEGFA -> KDR
# - HBC/PAMM1 -> trophoblast/core: IL1B/chemokines -> IL1R1/CCR/CXCR axes
#
# EXPECTED RESULT
# ---------------
# Targeted stromal -> trophoblast / endothelium and immune -> trophoblast
# interactions should be enriched in infected 24 hpi samples relative to UI.
# ============================================================

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("15_targeted_cellchat_axes starting.", log_file = log_file)

if (!requireNamespace("CellChat", quietly = TRUE)) {
  log_msg("CellChat not installed; skipping 15_targeted_cellchat_axes.", log_file = log_file)
  message("CellChat not installed; script skipped.")
} else {
  
  suppressPackageStartupMessages({
    library(CellChat)
  })
  
  fig_dir <- CFG$dirs$figures
  obj_dir <- CFG$dirs$objects
  tbl_dir <- CFG$dirs$tables
  
  # Prefer mapped object if available
  if (file.exists(file.path(obj_dir, "seu_with_architecture_transfer.qs"))) {
    seu <- qs::qread(file.path(obj_dir, "seu_with_architecture_transfer.qs"))
  } else {
    seu <- qs::qread(file.path(obj_dir, "seu_with_scores.qs"))
  }
  
  DefaultAssay(seu) <- "RNA"
  
  ct_col <- CFG$cols$cell_type
  inf_col <- CFG$cols$infection
  hpi_col <- CFG$cols$hpi
  
  # Build condition if absent
  if (!"condition" %in% colnames(seu@meta.data) &&
      all(c(inf_col, hpi_col) %in% colnames(seu@meta.data))) {
    seu$condition <- paste0(seu[[inf_col]][,1], "_", seu[[hpi_col]][,1])
  }
  
  # Use transferred labels when available
  if ("arch_predicted_label" %in% colnames(seu@meta.data)) {
    seu$comm_label <- ifelse(
      !is.na(seu$arch_predicted_label) & seu$arch_prediction_score_max >= 0.4,
      seu$arch_predicted_label,
      as.character(seu@meta.data[[ct_col]])
    )
  } else {
    seu$comm_label <- as.character(seu@meta.data[[ct_col]])
  }
  
  # ----------------------------------------------------------
  # Programmatic filtering: keep UI + 24 hpi only
  # ----------------------------------------------------------
  md <- seu@meta.data
  
  inf_vec <- as.character(md[[inf_col]])
  hpi_vec <- as.character(md[[hpi_col]])
  
  is_ui <- grepl("^UI$|uninfected|control", inf_vec, ignore.case = TRUE)
  is_24 <- hpi_vec %in% c("24", "24h", "24_h", "24 hpi")
  
  cells_keep <- rownames(md)[is_ui | is_24]
  
  if (length(cells_keep) == 0) {
    stop("No cells matched UI or 24 hpi for CellChat analysis.")
  }
  
  seu <- subset(seu, cells = cells_keep)
  
  # Optional focus on labels relevant to the targeted axes
  keep_labels <- unique(seu$comm_label)[
    grepl("FIB2|FIB|F$|F_|PV|Endo|EVT|VCT|HBC|PAMM", unique(seu$comm_label), ignore.case = TRUE)
  ]
  
  if (length(keep_labels) > 0) {
    cells_keep2 <- rownames(seu@meta.data)[seu$comm_label %in% keep_labels]
    seu <- subset(seu, cells = cells_keep2)
  }
  
  # CellChat works on normalized data; normalize if needed
  seu <- NormalizeData(seu, verbose = FALSE)
  
  # ----------------------------------------------------------
  # Build CellChat object
  # ----------------------------------------------------------
  cellchat <- createCellChat(
    object = seu,
    group.by = "comm_label",
    assay = "RNA"
  )
  
  cellchat@DB <- CellChatDB.human
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  qs::qsave(cellchat, file.path(obj_dir, "cellchat_targeted_24h.qs"), preset = "high")
  
  # ----------------------------------------------------------
  # Extract targeted interactions
  # ----------------------------------------------------------
  comm_all <- subsetCommunication(cellchat)
  
  target_ligands <- c("AREG", "HGF", "VEGFA", "IL1B", "CXCL8", "CCL20", "LGALS3")
  target_receptors <- c("EGFR", "MET", "KDR", "IL1R1", "CXCR1", "CXCR2", "CCR1", "CCR5")
  
  comm_target <- comm_all %>%
    filter(
      ligand %in% target_ligands |
        receptor %in% target_receptors |
        interaction_name %in% paste0(target_ligands, "_", target_receptors)
    )
  
  write.csv(comm_all, file.path(tbl_dir, "cellchat_all_24h_interactions.csv"), row.names = FALSE)
  write.csv(comm_target, file.path(tbl_dir, "cellchat_targeted_24h_interactions.csv"), row.names = FALSE)
  
  # ----------------------------------------------------------
  # Simple manuscript-ready plot
  # ----------------------------------------------------------
  if (nrow(comm_target) > 0) {
    plot_df <- comm_target %>%
      mutate(source_target = paste(source, "->", target)) %>%
      group_by(source_target, pathway_name) %>%
      summarise(prob = mean(prob, na.rm = TRUE), .groups = "drop")
    
    p <- ggplot(plot_df, aes(x = pathway_name, y = source_target, fill = prob)) +
      geom_tile(color = "white") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "Targeted CellChat axes at UI + 24 hpi",
        x = "Pathway",
        y = "Source -> Target",
        fill = "Mean\nprobability"
      )
    
    save_plot(file.path(fig_dir, "Fig36_Targeted_CellChat_axes.pdf"), p, w = 14, h = 8)
  } else {
    log_msg("No targeted CellChat interactions found for plotting.", log_file = log_file)
  }
  
  write_legend(
    "Fig36", "Targeted CellChat axes",
    hypothesis = "The Fn hypothesis is strengthened if stromal, vascular, and macrophage communication axes relevant to EA support or vascular stress are detectable in the infection object.",
    methods = "We ran CellChat on UI plus 24 hpi cells using original or transferred communication labels and extracted targeted ligand-receptor axes centered on AREG/HGF, VEGFA, IL1B, chemokines, and LGALS3-related signaling.",
    readout = "Evidence is strongest if stromal-to-trophoblast/endothelial and macrophage-to-core/trophoblast pathways are present with biologically plausible directionality.",
    interpretation_template = "- FIB/FIB2-like -> trophoblast via AREG/HGF supports a trophoblast-support niche.\n- FIB/FIB2-like -> endothelium via VEGFA -> KDR supports stromal-vascular coupling.\n- HBC/PAMM1 inflammatory axes support immune amplification of tissue injury.",
    outfile = file.path(CFG$dirs$legends, "Fig36_legend.md")
  )
  
  log_msg("15_targeted_cellchat_axes done.", log_file = log_file)
}