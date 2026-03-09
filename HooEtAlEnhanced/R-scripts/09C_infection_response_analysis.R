# ======================================================================
# 09C_infection_response_analysis.R
# Baseline infection response capacity analysis
#
# Analyzes the in vivo placenta's baseline capacity to respond to
# infection, based on pre-infection expression of:
# 1. Inflammatory cytokines/chemokines
# 2. Interferon-stimulated genes
# 3. TLR/PRR sensing machinery
# 4. Antimicrobial effectors
#
# Connects to Hoo et al. 2024 findings on placental infection response
# ======================================================================

source("config/config.R")
source("scripts/R/utils.R")

logfile <- file.path(DIR_LOGS, "09C_infection_response.log")
log_msg("=== 09C: Infection Response Analysis ===", logfile)

# --- Infection response gene sets ---
INFECTION_GENESETS <- list(
  # From Hoo et al. 2024 - genes upregulated upon infection
  General_Inflammatory = c("CXCL3", "CXCL8", "CCL20", "CCL4", "CCL3", 
                           "IL1B", "TNF", "NFKBIA", "TNFAIP3"),
  IFN_Response = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", 
                   "OAS1", "OAS2", "OAS3", "IFI6", "IFI27"),
  HBC_Antimicrobial = c("IL1B", "TNF", "CXCL8", "CCL3", "CCL4",
                        "GBP1", "GBP2", "GBP4", "GBP5", "IDO1"),
  TLR_Sensing = c("TLR2", "TLR4", "TLR5", "TLR7", "TLR9",
                  "MYD88", "TICAM1", "IRF3", "IRF7"),
  Complement = c("C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B", "CFB"),
  Phagocytosis = c("CD68", "CD163", "FCGR1A", "FCGR2A", "FCGR3A",
                   "MSR1", "MARCO", "MRC1"),
  
  # Fn-specific response genes (predicted)
  Fn_Predicted_Response = c("IL1B", "TNF", "CXCL8", "CCL2", "CCL3",
                            "NFKBIA", "TNFAIP3", "PTGS2", "IL6",
                            "MMP2", "MMP9", "TIMP1")
)

# --- Load data ---
st <- tryCatch(readRDS(file.path(DIR_OBJS, "slidetags_harmonized.rds")),
               error = function(e) readRDS(PATH_SLIDETAGS_RDS))

st <- safe_join_layers(st)

# --- Score all infection response modules ---
for (gs_name in names(INFECTION_GENESETS)) {
  genes <- INFECTION_GENESETS[[gs_name]]
  score_name <- paste0("inf_", tolower(gsub("[^A-Za-z0-9]", "_", gs_name)))
  st <- add_module_score_safe(st, genes, score_name)
  log_msg(paste("Scored:", gs_name), logfile)
}

# --- Compute immune readiness index ---
md <- st@meta.data
inf_cols <- grep("^inf_", colnames(md), value = TRUE)

if (length(inf_cols) >= 3) {
  # Normalize each score to 0-1
  for (col in inf_cols) {
    vals <- md[[col]]
    if (!all(is.na(vals))) {
      rng <- range(vals, na.rm = TRUE)
      if (diff(rng) > 0) {
        md[[paste0(col, "_norm")]] <- (vals - rng[1]) / diff(rng)
      }
    }
  }
  
  # Composite readiness
  norm_cols <- grep("_norm$", colnames(md), value = TRUE)
  norm_cols <- norm_cols[grepl("^inf_", norm_cols)]
  
  if (length(norm_cols) > 0) {
    md$immune_readiness_index <- rowMeans(md[, norm_cols, drop = FALSE], na.rm = TRUE)
    log_msg(paste("Immune readiness range:", 
                  round(min(md$immune_readiness_index, na.rm=TRUE), 3), "to",
                  round(max(md$immune_readiness_index, na.rm=TRUE), 3)), logfile)
  }
  
  st@meta.data <- md
}

# --- Summary by cell type ---
celltype_col <- detect_celltype_column(st)

if (!is.null(celltype_col)) {
  score_cols <- c(inf_cols, "immune_readiness_index")
  score_cols <- score_cols[score_cols %in% colnames(st@meta.data)]
  
  readiness_summary <- st@meta.data %>%
    dplyr::group_by(!!sym(celltype_col)) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      dplyr::across(dplyr::all_of(score_cols), 
                    ~mean(.x, na.rm = TRUE), 
                    .names = "mean_{.col}"),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(mean_immune_readiness_index))
  
  write.csv(readiness_summary,
            file.path(DIR_TABLES, "infection_readiness_by_celltype.csv"),
            row.names = FALSE)
  
  cat("\n=== Immune Readiness Ranking ===\n")
  print(readiness_summary[, c(celltype_col, "n_cells", "mean_immune_readiness_index")])
}

# --- Save ---
saveRDS(st, file.path(DIR_OBJS, "slidetags_infection_scored.rds"))

log_msg("=== 09C complete ===", logfile)
cat("\n=== 09C: Infection Response Analysis COMPLETE ===\n")
