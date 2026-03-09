source("R/config.R")
source("R/helpers_io.R")

ensure_dir(CFG$dirs$data)
ensure_dir(CFG$dirs$outputs)
ensure_dir(CFG$dirs$objects)
log_file <- file.path(CFG$dirs$logs, "pipeline.log")

log_msg("01_load_make_seurat starting.", log_file = log_file)

seu <- NULL

if (file.exists(CFG$data$seu_qs_path)) {
  log_msg("Loading Seurat from qs: ", CFG$data$seu_qs_path, log_file = log_file)
  seu <- qs::qread(CFG$data$seu_qs_path)
} else if (file.exists(CFG$data$seu_rds_path)) {
  log_msg("Loading Seurat from rds: ", CFG$data$seu_rds_path, log_file = log_file)
  seu <- readRDS(CFG$data$seu_rds_path)
}

if (is.null(seu)) {
  if (is.null(CFG$data$h5ad_path)) {
    stop("No Seurat object found and CFG$data$h5ad_path is NULL. Put seu.rds/seu.qs in data/ or set h5ad_path.")
  }
  if (!requireNamespace("zellkonverter", quietly = TRUE)) {
    stop("zellkonverter not installed. Install Bioconductor package zellkonverter.")
  }
  log_msg("Reading h5ad: ", CFG$data$h5ad_path, log_file = log_file)
  sce <- zellkonverter::readH5AD(CFG$data$h5ad_path)
  seu <- Seurat::as.Seurat(sce, counts = "counts", data = "logcounts")
}

stopifnot(inherits(seu, "Seurat"))
DefaultAssay(seu) <- "RNA"

needed <- unlist(CFG$cols)
missing <- setdiff(needed, colnames(seu@meta.data))
if (length(missing)) {
  log_msg("WARNING: missing metadata columns: ", paste(missing, collapse=", "), log_file = log_file)
}

seu$hpi_num <- readr::parse_number(as.character(seu[[CFG$cols$hpi]][,1]))
seu$condition <- paste0(seu[[CFG$cols$infection]][,1], "_", seu[[CFG$cols$hpi]][,1])

saveRDS(seu, file = file.path(CFG$dirs$objects, "seu_clean.rds"))
qs::qsave(seu, file = file.path(CFG$dirs$objects, "seu_clean.qs"))

log_msg("Saved objects to outputs/objects. 01_load_make_seurat done.", log_file = log_file)
