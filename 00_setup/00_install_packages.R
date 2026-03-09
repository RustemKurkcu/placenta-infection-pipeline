# 00_setup/00_install_packages.R
# Run once in a fresh R environment

pkgs <- c(
  "Seurat", "SeuratObject", "Matrix", "dplyr", "ggplot2", "patchwork",
  "readr", "stringr", "tibble", "data.table", "qs",
  "cowplot", "RColorBrewer"
)

to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)

message("Done. Restart R, then run scripts/01_load_make_seurat.R")
