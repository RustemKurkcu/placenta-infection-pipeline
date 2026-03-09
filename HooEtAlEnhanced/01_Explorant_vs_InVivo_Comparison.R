# ============================================================================
# Explant vs In Vivo Placenta Comparison
# ============================================================================
# Comprehensive R script for comparing explant (Hoo et al. 2024) and 
# in vivo (Greenbaum et al.) placental data
#
# Focus: Uninfected cells, F. nucleatum vulnerability, spatial behavior
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(patchwork)
  library(scales)
})

# ============================================================================
# CONFIGURATION
# ============================================================================

CONFIG <- list(
  # Input paths
  explant_metadata = "/workspace/metadata.csv",
  invivo_celltype_props = "hPlacenta-architecture/output/tables/celltype_proportions_by_week_and_version.csv",
  
  # Output paths
  output_dir = "/workspace/analysis/R_scripts/output/",
  fig_dir = "/workspace/analysis/R_scripts/figures/",
  
  # Quality thresholds
  min_features = 500,
  max_mt_percent = 20
)

# Create output directories
dir.create(CONFIG$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(CONFIG$fig_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading data...\n")

# Load explant data
explant_meta <- read_csv(CONFIG$explant_metadata, show_col_types = FALSE)
explant_ui <- explant_meta %>% filter(infection == "UI")
cat(sprintf("  Explant: %s total cells, %s uninfected cells\n", 
            nrow(explant_meta), nrow(explant_ui)))

# Load in vivo data
invivo_props <- read_csv(CONFIG$invivo_celltype_props, show_col_types = FALSE)
invivo_agg <- invivo_props %>%
  group_by(celltype) %>%
  summarize(
    invivo_n = sum(n_cells),
    invivo_mean_prop = mean(prop),
    .groups = "drop"
  ) %>%
  arrange(desc(invivo_n))
cat(sprintf("  In Vivo: %s cells, %s cell types\n", 
            sum(invivo_agg$invivo_n), nrow(invivo_agg)))

# ============================================================================
# CELL TYPE HARMONIZATION
# ============================================================================

HARMONIZATION_MAP <- list(
  "HBC" = "Hofbauer cells",
  "HBC_p" = "Hofbauer cells",
  "F" = "FIB1",
  "F_p" = "FIB1",
  "F_sm" = "FIB2",
  "VCT" = "vCTB",
  "VCT_CCC" = "vCTB",
  "VCT_p" = "vCTB",
  "VCT_fusing" = "vCTB",
  "EVT_1" = "EVT1",
  "EVT_2" = "EVT",
  "iEVT" = "EVT-progenitor",
  "Endo_f" = "Endothelial",
  "PAMM1" = "STB-progenitor",
  "PV" = "STB"
)

# Apply harmonization
explant_counts <- explant_ui %>%
  count(cell_type, name = "explant_n") %>%
  mutate(
    explant_pct = (explant_n / sum(explant_n)) * 100,
    harmonized_ct = unlist(HARMONIZATION_MAP[cell_type])
  )

# Merge datasets
comparison <- invivo_agg %>%
  left_join(explant_counts, by = c("celltype" = "harmonized_ct")) %>%
  mutate(
    explant_n = ifelse(is.na(explant_n), 0, explant_n),
    explant_pct = ifelse(is.na(explant_pct), 0, explant_pct),
    ratio_explant_to_invivo = ifelse(
      invivo_mean_prop > 0,
      explant_pct / (invivo_mean_prop * 100),
      NA_real_
    )
  )

write_csv(comparison, file.path(CONFIG$output_dir, "explant_vs_invivo_comparison.csv"))

# ============================================================================
# FIGURE 1: CELL TYPE COMPOSITION
# ============================================================================

cat("Creating Figure 1: Cell Type Composition...\n")

fig1_a <- ggplot(explant_counts, aes(x = reorder(cell_type, -explant_n), y = explant_n)) +
  geom_bar(stat = "identity", fill = "#3498db") +
  coord_flip() +
  labs(
    title = "A. Explant Model (Uninfected)",
    subtitle = "Hoo et al. 2024",
    x = "Cell Type",
    y = "Number of Cells"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_text(size = 9))

fig1_b <- ggplot(invivo_agg %>% slice_max(invivo_n, n = 10), 
                 aes(x = reorder(celltype, -invivo_n), y = invivo_n)) +
  geom_bar(stat = "identity", fill = "#2ecc71") +
  coord_flip() +
  labs(
    title = "B. In Vivo Placenta",
    subtitle = "Greenbaum et al.",
    x = "Cell Type",
    y = "Number of Cells"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_text(size = 9))

fig1 <- fig1_a + fig1_b + plot_layout(ncol = 2)
ggsave(file.path(CONFIG$fig_dir, "Fig01_CellType_Composition.pdf"), 
       fig1, width = 14, height = 6)

# ============================================================================
# FIGURE 2: FN VULNERABILITY SCORES
# ============================================================================

cat("Creating Figure 2: F. nucleatum Vulnerability...\n")

# Calculate vulnerability scores
score_summary <- explant_ui %>%
  group_by(cell_type) %>%
  summarize(
    n_cells = n(),
    glyco_mean = mean(GlycoScore, na.rm = TRUE),
    adhesion_mean = mean(AdhesionScore, na.rm = TRUE),
    innate_mean = mean(InnateScore, na.rm = TRUE),
    cytotoxic_mean = mean(CytotoxicScore, na.rm = TRUE),
    fn_vulnerability = (
      0.3 * mean(GlycoScore, na.rm = TRUE) +
      0.3 * mean(AdhesionScore, na.rm = TRUE) +
      0.4 * (1 - mean(InnateScore, na.rm = TRUE))
    ),
    .groups = "drop"
  ) %>%
  arrange(fn_vulnerability)

fig2_a <- ggplot(score_summary, 
                 aes(x = reorder(cell_type, glyco_mean), y = glyco_mean)) +
  geom_bar(stat = "identity", aes(fill = glyco_mean)) +
  scale_fill_gradient(low = "#3498db", high = "#e74c3c") +
  coord_flip() +
  labs(
    title = "A. Glycosylation Score",
    subtitle = "Fn Adhesion Potential",
    x = "Cell Type",
    y = "Mean GlycoScore"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

fig2_b <- ggplot(score_summary, 
                 aes(x = reorder(cell_type, adhesion_mean), y = adhesion_mean)) +
  geom_bar(stat = "identity", aes(fill = adhesion_mean)) +
  scale_fill_gradient(low = "#3498db", high = "#e74c3c") +
  coord_flip() +
  labs(
    title = "B. Adhesion Score",
    subtitle = "Fn Binding Potential",
    x = "Cell Type",
    y = "Mean AdhesionScore"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

fig2_c <- ggplot(score_summary, 
                 aes(x = reorder(cell_type, innate_mean), y = innate_mean)) +
  geom_bar(stat = "identity", aes(fill = innate_mean)) +
  scale_fill_gradient(low = "#e74c3c", high = "#2ecc71") +
  coord_flip() +
  labs(
    title = "C. Innate Immune Score",
    subtitle = "Defense Response",
    x = "Cell Type",
    y = "Mean InnateScore"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

fig2_d <- ggplot(score_summary, 
                 aes(x = reorder(cell_type, fn_vulnerability), y = fn_vulnerability)) +
  geom_bar(stat = "identity", aes(fill = fn_vulnerability)) +
  scale_fill_gradient(low = "#3498db", high = "#e74c3c") +
  coord_flip() +
  labs(
    title = "D. Combined Fn Vulnerability Index",
    subtitle = "0.3*Glyco + 0.3*Adhesion + 0.4*(1-Innate)",
    x = "Cell Type",
    y = "Fn Vulnerability Index"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

fig2 <- (fig2_a + fig2_b) / (fig2_c + fig2_d)
ggsave(file.path(CONFIG$fig_dir, "Fig02_FN_Vulnerability.pdf"), 
       fig2, width = 12, height = 12)

# ============================================================================
# FIGURE 3: MODEL VALIDATION
# ============================================================================

cat("Creating Figure 3: Model Validation...\n")

# Quality metrics summary
qc_summary <- data.frame(
  Metric = c("Median Genes", "Median UMIs (K)", "MT% (<5%)", 
             "Cell Type Diversity", "Donor Count"),
  Value = c(
    median(explant_ui$nFeature_RNA),
    median(explant_ui$nCount_RNA) / 1000,
    100 - (mean(explant_ui$percent.mt) * 100),
    length(unique(explant_ui$cell_type)),
    length(unique(explant_ui$donor_id))
  )
)

fig3_a <- ggplot(qc_summary, aes(x = reorder(Metric, -Value), y = Value)) +
  geom_bar(stat = "identity", fill = "#9b59b6") +
  coord_flip() +
  labs(
    title = "A. Quality Metrics",
    x = "Metric",
    y = "Value"
  ) +
  theme_minimal()

# Proportion comparison scatter
matched_comparison <- comparison %>%
  filter(explant_n > 0) %>%
  filter(!is.na(ratio_explant_to_invivo))

fig3_b <- ggplot(matched_comparison, 
                 aes(x = invivo_mean_prop * 100, y = explant_pct)) +
  geom_point(size = 4, alpha = 0.7, color = "#3498db") +
  geom_text(aes(label = celltype), size = 3, vjust = -0.5, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    title = "B. Cell Type Proportion Comparison",
    x = "In Vivo Proportion (%)",
    y = "Explant Proportion (%)"
  ) +
  theme_minimal() +
  xlim(0, 30) +
  ylim(0, 30)

fig3 <- fig3_a + fig3_b + plot_layout(ncol = 2, widths = c(1, 2))
ggsave(file.path(CONFIG$fig_dir, "Fig03_Model_Validation.pdf"), 
       fig3, width = 14, height = 6)

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("Saving results...\n")

# Save score summary
write_csv(score_summary, file.path(CONFIG$output_dir, "fn_vulnerability_summary.csv"))

# Create summary report
summary_report <- list(
  dataset = list(
    explant_total = nrow(explant_meta),
    explant_uninfected = nrow(explant_ui),
    explant_celltypes = length(unique(explant_ui$cell_type)),
    invivo_total = sum(invivo_agg$invivo_n),
    invivo_celltypes = nrow(invivo_agg),
    matched_celltypes = length(unique(comparison$celltype[comparison$explant_n > 0]))
  ),
  quality = list(
    median_genes = median(explant_ui$nFeature_RNA),
    median_umis = median(explant_ui$nCount_RNA),
    median_mt = median(explant_ui$percent.mt) * 100
  ),
  vulnerability = list(
    most_vulnerable = score_summary$cell_type[which.max(score_summary$fn_vulnerability)],
    least_vulnerable = score_summary$cell_type[which.min(score_summary$fn_vulnerability)]
  )
)

saveRDS(summary_report, file.path(CONFIG$output_dir, "summary_report.rds"))

cat("\nAnalysis complete!\n")
cat(sprintf("Figures saved to: %s\n", CONFIG$fig_dir))
cat(sprintf("Results saved to: %s\n", CONFIG$output_dir))