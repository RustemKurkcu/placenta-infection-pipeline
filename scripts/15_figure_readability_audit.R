source("R/config.R")
source("R/helpers_io.R")

suppressPackageStartupMessages({
  library(dplyr)
})

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("15_figure_readability_audit starting.", log_file = log_file)

fig_dir <- CFG$dirs$figures
leg_dir <- CFG$dirs$legends
tbl_dir <- CFG$dirs$tables

figs <- list.files(fig_dir, pattern = "\\.(pdf|png)$", full.names = TRUE)
if (length(figs) == 0) {
  log_msg("No figures found in ", fig_dir, ".", log_file = log_file)
  write.csv(data.frame(note = "No figures available for audit"), file.path(tbl_dir, "figure_readability_audit.csv"), row.names = FALSE)
  quit(save = "no")
}

info <- file.info(figs)
audit <- data.frame(
  file = basename(figs),
  path = figs,
  size_bytes = info$size,
  ext = tools::file_ext(figs),
  stringsAsFactors = FALSE
)

audit <- audit %>%
  mutate(
    size_kb = round(size_bytes / 1024, 2),
    readability_flag = case_when(
      ext == "pdf" & size_kb < 20 ~ "likely_low_detail",
      ext == "png" & size_kb < 80 ~ "likely_low_detail",
      TRUE ~ "ok"
    ),
    has_matching_legend = file.exists(file.path(leg_dir, gsub("\\.(pdf|png)$", "_legend.md", file)))
  ) %>%
  arrange(readability_flag, size_kb)

write.csv(audit, file.path(tbl_dir, "figure_readability_audit.csv"), row.names = FALSE)

flagged <- audit %>% filter(readability_flag != "ok" | !has_matching_legend)
write.csv(flagged, file.path(tbl_dir, "figure_readability_flagged.csv"), row.names = FALSE)

md <- c(
  "# Figure readability audit",
  "",
  paste0("Total figures audited: ", nrow(audit)),
  paste0("Flagged figures: ", nrow(flagged)),
  "",
  "## Flag rules",
  "- PDF < 20KB => likely low visual detail/possibly unreadable.",
  "- PNG < 80KB => likely low visual detail/possibly unreadable.",
  "- Missing legend file also flagged.",
  ""
)

if (nrow(flagged) > 0) {
  md <- c(md, "## Flagged figures", "")
  for (i in seq_len(nrow(flagged))) {
    md <- c(md, paste0("- ", flagged$file[i], " | ", flagged$readability_flag[i], " | size_kb=", flagged$size_kb[i], " | legend=", flagged$has_matching_legend[i]))
  }
} else {
  md <- c(md, "No flagged figures by heuristic rules.")
}

writeLines(md, con = file.path(tbl_dir, "figure_readability_audit.md"))
log_msg("15_figure_readability_audit done.", log_file = log_file)
