source("R/config.R")
source("R/helpers_io.R")

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("05_susceptibility_severity_models starting.", log_file = log_file)

seu <- qs::qread(file.path(CFG$dirs$objects, "seu_with_scores.qs"))

# Define sample_id (edit fields if needed)
seu$sample_id <- paste(
  seu[[CFG$cols$donor_id]][,1],
  seu[[CFG$cols$infection]][,1],
  seu[[CFG$cols$hpi]][,1],
  seu[[CFG$cols$stage]][,1],
  sep="|"
)

df_sample <- seu@meta.data %>%
  group_by(sample_id,
           donor_id = .data[[CFG$cols$donor_id]],
           infection = .data[[CFG$cols$infection]],
           hpi = .data[[CFG$cols$hpi]],
           stage = .data[[CFG$cols$stage]]) %>%
  summarise(
    n_cells = n(),
    glyco    = mean(GlycoScore, na.rm = TRUE),
    adhesion = mean(AdhesionScore, na.rm = TRUE),
    innate   = mean(InnateScore, na.rm = TRUE),
    .groups = "drop"
  )

ensure_dir(CFG$dirs$tables)
write.csv(df_sample, file.path(CFG$dirs$tables, "sample_level_scores.csv"), row.names = FALSE)

# Correlations within infection x time

df_cor <- df_sample %>%
  group_by(infection, hpi) %>%
  summarise(
    cor_glyco_innate = suppressWarnings(cor(glyco, innate, use="complete.obs")),
    cor_adh_innate   = suppressWarnings(cor(adhesion, innate, use="complete.obs")),
    n_samples = n(),
    .groups = "drop"
  )

write.csv(df_cor, file.path(CFG$dirs$tables, "cor_susceptibility_vs_severity.csv"), row.names = FALSE)

write_legend(
  "Fig12", "Susceptibility vs severity at sample level (correlations)",
  hypothesis = "Baseline/average entry-susceptibility (glyco/adhesion) predicts inflammatory severity (innate score), possibly in a pathogen- and time-dependent manner.",
  methods = "Aggregate per-sample mean scores, then compute correlations within infection × time strata.",
  readout = "Positive correlation suggests susceptibility→severity coupling; negative/zero suggests decoupled programs or confounding.",
  interpretation_template = "- If coupling exists for bacterial infection (Lm), it may be informative for Fn-like hypotheses.\n- If absent, consider cell-type specific pseudo-bulk (trophoblast-only susceptibility; macrophage-only severity).",
  outfile = file.path(CFG$dirs$legends, "Fig12_legend.md")
)

log_msg("05_susceptibility_severity_models done.", log_file = log_file)
