source("R/config.R")
source("R/helpers_io.R")

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
})

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("09_reproducibility_report starting.", log_file = log_file)

obj_path <- file.path(CFG$dirs$objects, "seu_with_scores.qs")
if (!file.exists(obj_path)) stop("Missing object: ", obj_path)

seu <- qs::qread(obj_path)
md <- seu@meta.data

ensure_dir(CFG$dirs$tables)
ensure_dir(CFG$dirs$logs)

# Replicate table by infection/time/donor
rep_tab <- md %>%
  count(infection = .data[[CFG$cols$infection]],
        hpi = .data[[CFG$cols$hpi]],
        donor_id = .data[[CFG$cols$donor_id]],
        name = "n_cells") %>%
  arrange(infection, hpi, donor_id)

write.csv(rep_tab, file.path(CFG$dirs$tables, "replicate_cell_counts_by_condition.csv"), row.names = FALSE)

# Condition-level donor counts
donor_tab <- rep_tab %>%
  group_by(infection, hpi) %>%
  summarise(n_donors = n_distinct(donor_id), total_cells = sum(n_cells), .groups = "drop")

write.csv(donor_tab, file.path(CFG$dirs$tables, "donor_counts_by_condition.csv"), row.names = FALSE)

# Session info for reproducibility
sink(file.path(CFG$dirs$logs, "session_info.txt"))
print(sessionInfo())
sink()

log_msg("09_reproducibility_report done.", log_file = log_file)
