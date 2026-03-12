source("R/config.R")
source("R/helpers_io.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(qs)
})

# ============================================================
# 12_pathogen_specificity_matrix.R
#
# HYPOTHESIS
# ----------
# Not all placental pathogens or oral bacteria should be expected to drive
# the same PE-like vascular syndrome. This script creates a transparent,
# manuscript-ready comparison matrix that separates shared features from
# hypothesized Fn-specific mechanisms.
# ============================================================

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("12_pathogen_specificity_matrix starting.", log_file = log_file)

fig_dir <- CFG$dirs$figures
tbl_dir <- CFG$dirs$tables
obj_dir <- CFG$dirs$objects

mat <- tibble::tribble(
  ~pathogen, ~placental_crossing, ~ethanolamine_utilization, ~BMC_support, ~H2S_generation, ~endothelial_dysfunction_evidence, ~antiangiogenic_FLT1_link, ~preeclampsia_association, ~notes,
  "Fusobacterium_nucleatum", 1, 1, 1, 1, 1, 0.5, 1, "Central hypothesis organism; strongest MegL/H2S + EA/BMC logic, but direct PE causality is not yet proven.",
  "Listeria_monocytogenes", 1, 1, 1, 0, 0.5, 0, 0, "Best EA-positive / MegL-negative comparator.",
  "Porphyromonas_gingivalis", 0.5, 0, 0, 1, 1, 0.5, 1, "Best oral/endothelial comparator for host dysfunction.",
  "Toxoplasma_gondii", 1, 0, 0, 0, 0, 0, 0, "Useful negative control for broad inflammation vs vascular/EA mechanisms.",
  "Plasmodium_falciparum", 0, 0, 0, 0, 1, 1, 0.5, "Positive control for placental vascular stress and FLT1/PGF disruption."
)

write.csv(mat, file.path(tbl_dir, "pathogen_specificity_matrix.csv"), row.names = FALSE)
qs::qsave(mat, file.path(obj_dir, "pathogen_specificity_matrix.qs"), preset = "high")

long_df <- mat %>% select(-notes) %>% pivot_longer(-pathogen, names_to = "feature", values_to = "score")
p <- ggplot(long_df, aes(x = feature, y = pathogen, fill = score)) +
  geom_tile(color = "white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Pathogen specificity matrix for the Fn-PE hypothesis", x = NULL, y = NULL, fill = "evidence\nscore")

save_plot(file.path(fig_dir, "Fig34_Pathogen_specificity_matrix.pdf"), p, w = 14, h = 5)

write_legend(
  "Fig34", "Pathogen specificity matrix",
  hypothesis = "Fn should differ from other pathogens by combining a placental niche, EA/BMC support, and a putative H2S-linked vascular toxicity program.",
  methods = "We manually curated a comparative matrix across pathogens used or discussed in the project and scored each one for placental crossing, EA utilization, BMC support, H2S generation, endothelial dysfunction evidence, FLT1-linked evidence, and PE association.",
  readout = "The matrix should clarify which features are shared and which are specific to the proposed Fn MegL toxic-switch model.",
  interpretation_template = "- Listeria positive in EA/BMC but negative in MegL-like H2S supports specificity.\n- P. falciparum positive in vascular stress but negative in EA biology supports pathway decomposition.\n- Oral/endothelial comparators such as P. gingivalis strengthen the host-damage axis without proving the EA switch.",
  outfile = file.path(CFG$dirs$legends, "Fig34_legend.md")
)

log_msg("12_pathogen_specificity_matrix done.", log_file = log_file)
