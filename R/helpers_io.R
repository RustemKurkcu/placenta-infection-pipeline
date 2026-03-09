ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_msg <- function(..., log_file = NULL) {
  msg <- paste0("[", timestamp(), "] ", paste(..., collapse = ""))
  message(msg)
  if (!is.null(log_file)) {
    ensure_dir(dirname(log_file))
    cat(msg, "\n", file = log_file, append = TRUE)
  }
}

save_plot <- function(path, p, w = 10, h = 7, dpi = 300) {
  ensure_dir(dirname(path))
  ggsave(filename = path, plot = p, width = w, height = h, units = "in", dpi = dpi, limitsize = FALSE)
  invisible(path)
}

write_legend <- function(fig_id, title, hypothesis, methods, readout, interpretation_template,
                         outfile = NULL) {
  if (is.null(outfile)) outfile <- file.path(CFG$dirs$legends, paste0(fig_id, "_legend.md"))
  ensure_dir(dirname(outfile))
  txt <- paste(
    paste0("# ", fig_id, ": ", title),
    "",
    "## Hypothesis / idea being tested",
    hypothesis,
    "",
    "## Methods (what we did)",
    methods,
    "",
    "## Readout (what to look for)",
    readout,
    "",
    "## Interpretation template",
    interpretation_template,
    "",
    sep="\n"
  )
  writeLines(txt, con = outfile, useBytes = TRUE)
  invisible(outfile)
}
