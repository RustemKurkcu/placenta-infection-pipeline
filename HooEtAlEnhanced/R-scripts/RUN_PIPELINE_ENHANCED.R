# ======================================================================
# RUN_PIPELINE_ENHANCED.R
# Enhanced pipeline including Fn vulnerability and organoid comparison
# ======================================================================

source("config/config.R")

message("\n=== RUN_PIPELINE_ENHANCED ===")
message("Project root: ", getwd())
message("Seurat version: ", as.character(utils::packageVersion("Seurat")))

# --- Core Pipeline (Minimal) ---
message("\n--- Phase 1: Core Pipeline ---")
source("scripts/02_preprocess/02A_preprocess_multiome_reference.R")
source("scripts/03_mapping/03A_map_slidetags_to_multiome.R")
source("scripts/03_mapping/03B_map_starmap_to_multiome.R")
source("scripts/03_mapping/03C_harmonize_celltype_labels.R")

# --- Timecourse & Spatial ---
message("\n--- Phase 2: Timecourse & Spatial ---")
source("scripts/04_timecourse/04A_gene_of_interest_timecourse.R")
source("scripts/05_spatial/05A_spatial_overview_plots.R")
source("scripts/05_spatial/05B_neighborhood_enrichment.R")
source("scripts/05_spatial/05C_permissiveness_score_maps.R")

# --- Cell Communication ---
message("\n--- Phase 3: Cell Communication ---")
source("scripts/06_cell_communication/06B_simple_LR_scoring.R")

# --- Metagene Analysis ---
message("\n--- Phase 4: Metagene Analysis ---")
source("scripts/08_metagenes/08A_housekeeping_diagnostics.R")
source("scripts/08_metagenes/08B_metagene_module_discovery.R")

# --- NEW: Fn Vulnerability & Organoid Comparison ---
message("\n--- Phase 5: Fn Vulnerability & Organoid Analysis ---")
source("scripts/09_fn_analysis/09A_fn_vulnerability_scoring.R")
source("scripts/09_fn_analysis/09B_organoid_invivo_comparison.R")
source("scripts/09_fn_analysis/09C_infection_response_analysis.R")

# --- Export ---
message("\n--- Phase 6: Export ---")
source("scripts/07_export/07A_export_shareable_outputs.R")

message("\n=== DONE (ENHANCED) ===\n")
