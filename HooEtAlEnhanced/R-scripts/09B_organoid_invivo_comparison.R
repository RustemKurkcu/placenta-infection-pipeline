# ======================================================================
# 09B_organoid_invivo_comparison.R
# Compare organoid/explant gene expression with in vivo placenta
#
# This script provides a framework for comparing:
# 1. Cell type marker fidelity
# 2. Pathway conservation
# 3. EA metabolism landscape
# 4. Infection response capacity
#
# Requires: Hoo et al. 2024 data (E-MTAB-12795) or summary statistics
# ======================================================================

source("config/config.R")
source("scripts/R/utils.R")

logfile <- file.path(DIR_LOGS, "09B_organoid_comparison.log")
log_msg("=== 09B: Organoid vs In Vivo Comparison ===", logfile)

# --- Cell type mapping between datasets ---
CELLTYPE_MAP <- list(
  # In vivo (Greenbaum) -> Organoid/Explant (Hoo et al.)
  "vCTB" = c("VCT", "VCT_fusing", "VCT_p", "VCT_CCC"),
  "STB" = c("SCT"),
  "EVT" = c("iEVT", "EVT_1", "EVT_2"),
  "Hofbauer cells" = c("HBC", "HBC_p"),
  "FIB2" = c("F", "F_p"),
  "Endothelial" = c("Endo_f")
)

# --- Canonical markers for cross-platform validation ---
CANONICAL_MARKERS <- list(
  VCT = c("KRT7", "KRT8", "KRT18", "TACSTD2", "EPCAM", "TP63", "TEAD4", "CDH1"),
  SCT = c("CSH1", "CSH2", "CYP19A1", "PSG1", "CGA", "KISS1", "GCM1"),
  EVT = c("HLA-G", "ITGA1", "MMP2", "MMP9", "FN1", "ADAM12", "ITGA5"),
  HBC = c("CD68", "CD163", "CSF1R", "C1QA", "C1QB", "C1QC", "SPP1", "F13A1"),
  Fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "VIM"),
  Endothelial = c("PECAM1", "VWF", "KDR", "ESAM", "CDH5", "RAMP2")
)

# --- Load in vivo data ---
st <- tryCatch(readRDS(file.path(DIR_OBJS, "slidetags_harmonized.rds")),
               error = function(e) readRDS(PATH_SLIDETAGS_RDS))

log_msg(paste("Loaded in vivo data:", ncol(st), "cells"), logfile)

# --- Compute per-cell-type average expression ---
compute_celltype_profiles <- function(obj, celltype_col = NULL) {
  if (is.null(celltype_col)) {
    celltype_col <- detect_celltype_column(obj)
  }
  
  obj <- safe_join_layers(obj)
  
  # Get expression matrix
  expr_mat <- GetAssayData(obj, layer = "data")
  
  # Group by cell type
  cell_types <- obj@meta.data[[celltype_col]]
  unique_cts <- unique(cell_types)
  
  profiles <- list()
  for (ct in unique_cts) {
    cells <- which(cell_types == ct)
    if (length(cells) >= 5) {
      profiles[[ct]] <- list(
        mean_expr = Matrix::rowMeans(expr_mat[, cells, drop = FALSE]),
        detect_rate = Matrix::rowMeans(expr_mat[, cells, drop = FALSE] > 0),
        n_cells = length(cells)
      )
    }
  }
  
  return(profiles)
}

invivo_profiles <- compute_celltype_profiles(st)
log_msg(paste("Computed profiles for", length(invivo_profiles), "cell types"), logfile)

# --- Marker fidelity analysis ---
marker_fidelity <- function(profiles, markers_list) {
  results <- data.frame()
  
  for (ct_name in names(markers_list)) {
    markers <- markers_list[[ct_name]]
    
    for (gene in markers) {
      for (profile_name in names(profiles)) {
        prof <- profiles[[profile_name]]
        if (gene %in% names(prof$mean_expr)) {
          expr_val <- prof$mean_expr[gene]
          detect_val <- prof$detect_rate[gene]
          
          # Compute specificity vs other cell types
          other_means <- sapply(profiles[names(profiles) != profile_name], function(p) {
            if (gene %in% names(p$mean_expr)) p$mean_expr[gene] else 0
          })
          other_mean <- mean(other_means, na.rm = TRUE)
          specificity <- expr_val / (other_mean + 0.001)
          
          results <- rbind(results, data.frame(
            marker_set = ct_name,
            gene = gene,
            cell_type = profile_name,
            mean_expr = expr_val,
            detect_rate = detect_val,
            specificity = specificity,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  return(results)
}

fidelity_results <- marker_fidelity(invivo_profiles, CANONICAL_MARKERS)
write.csv(fidelity_results, 
          file.path(DIR_TABLES, "marker_fidelity_invivo.csv"),
          row.names = FALSE)
log_msg("Saved marker fidelity results", logfile)

# --- EA metabolism comparison ---
EA_GENES <- list(
  release = c("PLD1", "GDPD1", "GDPD5", "FAAH", "PLA2G6"),
  consumption = c("ETNK1", "ETNK2", "PCYT2", "SELENOI", "CEPT1"),
  synthesis = c("PCYT2", "SELENOI", "CEPT1", "CHPT1"),
  remodeling = c("LPCAT1", "LPCAT2", "LPCAT3", "LPCAT4")
)

ea_results <- data.frame()
for (ct_name in names(invivo_profiles)) {
  prof <- invivo_profiles[[ct_name]]
  for (pathway in names(EA_GENES)) {
    genes <- EA_GENES[[pathway]]
    found <- genes[genes %in% names(prof$mean_expr)]
    if (length(found) > 0) {
      total_expr <- sum(prof$mean_expr[found])
      ea_results <- rbind(ea_results, data.frame(
        cell_type = ct_name,
        pathway = pathway,
        n_genes = length(found),
        total_expr = total_expr,
        mean_expr = total_expr / length(found),
        stringsAsFactors = FALSE
      ))
    }
  }
}

write.csv(ea_results,
          file.path(DIR_TABLES, "ea_metabolism_by_celltype.csv"),
          row.names = FALSE)
log_msg("Saved EA metabolism results", logfile)

# --- Generate comparison plots ---
library(ggplot2)
library(tidyr)

# Marker fidelity heatmap
if (nrow(fidelity_results) > 0) {
  # Filter to target cell types
  target_pairs <- data.frame(
    marker_set = c("VCT", "SCT", "EVT", "HBC", "Fibroblast", "Endothelial"),
    target_ct = c("vCTB", "STB", "EVT", "Hofbauer cells", "FIB2", "Endothelial"),
    stringsAsFactors = FALSE
  )
  
  p_fidelity <- ggplot(fidelity_results, 
                        aes(x = cell_type, y = gene, fill = log1p(mean_expr))) +
    geom_tile(color = "white", linewidth = 0.3) +
    facet_wrap(~marker_set, scales = "free_y", ncol = 2) +
    scale_fill_viridis_c(option = "magma", name = "log(expr+1)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7),
          strip.text = element_text(face = "bold")) +
    labs(title = "Marker Gene Expression Across Cell Types (In Vivo)",
         subtitle = "Validating cell type identity for organoid comparison",
         x = "Cell Type", y = "Marker Gene")
  
  ggsave(file.path(DIR_FIGURES, "marker_fidelity_heatmap.png"),
         p_fidelity, width = 14, height = 10, dpi = 300)
  log_msg("Saved marker fidelity heatmap", logfile)
}

# EA metabolism bar plot
if (nrow(ea_results) > 0) {
  p_ea <- ggplot(ea_results, aes(x = reorder(cell_type, -mean_expr), 
                                  y = mean_expr, fill = pathway)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "EA Metabolism Pathway Expression by Cell Type",
         subtitle = "Release vs Consumption determines free EA for F. nucleatum",
         x = "Cell Type", y = "Mean Expression", fill = "EA Pathway")
  
  ggsave(file.path(DIR_FIGURES, "ea_metabolism_comparison.png"),
         p_ea, width = 12, height = 7, dpi = 300)
  log_msg("Saved EA metabolism plot", logfile)
}

log_msg("=== 09B complete ===", logfile)
cat("\n=== 09B: Organoid vs In Vivo Comparison COMPLETE ===\n")
