suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(qs)
})

CFG <- list(
  project = list(name = "placenta_infection_scRNA"),
  dirs = list(
    base = normalizePath(getwd(), winslash = "/", mustWork = FALSE),
    data = "data",
    outputs = "outputs",
    figures = file.path("outputs", "figures"),
    legends = file.path("outputs", "legends"),
    tables  = file.path("outputs", "tables"),
    objects = file.path("outputs", "objects"),
    logs    = file.path("outputs", "logs")
  ),
  cols = list(
    cell_type = "cell_type",
    infection = "infection",
    hpi       = "hpi",
    donor_id  = "donor_id",
    stage     = "stage"
  ),
  data = list(
    # set ONE of these
    seu_rds_path = file.path("data", "seu.rds"),
    seu_qs_path  = file.path("data", "seu.qs"),
    h5ad_path    = NULL   # e.g., "data/your_file.h5ad" if starting from h5ad
  ),
  plotting = list(
    dpi = 300,
    pt_size = 0.15
  )
)
