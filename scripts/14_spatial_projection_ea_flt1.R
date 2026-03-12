source("R/config.R")
source("R/helpers_io.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(qs)
  library(Matrix)
  library(patchwork)
})

# ============================================================
# 14_spatial_projection_ea_flt1.R
#
# HYPOTHESIS
# ----------
# If EA availability and vascular-stress programs are spatially structured,
# the architecture atlas should reveal discrete niches where these programs
# co-localize or interface with trophoblast / endothelial neighborhoods.
# ============================================================

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("14_spatial_projection_ea_flt1 starting.", log_file = log_file)

fig_dir <- CFG$dirs$figures
obj_dir <- CFG$dirs$objects
tbl_dir <- CFG$dirs$tables

# Build spatial reference object
meta <- read.csv("/mnt/data/metadata.csv", stringsAsFactors = FALSE)
meta <- meta[-1, , drop = FALSE]
coords <- read.csv("/mnt/data/humanplacenta_spatial.csv", stringsAsFactors = FALSE)
coords <- coords[-1, , drop = FALSE]
expr <- read.csv("/mnt/data/humanplacenta_expression.csv.gz", check.names = FALSE)

genes <- expr[[1]]
mat <- as.matrix(expr[, -1, drop = FALSE])
rownames(mat) <- genes
mat <- Matrix(mat, sparse = TRUE)

common_cells <- Reduce(intersect, list(colnames(mat), meta$NAME, coords$NAME))
mat <- mat[, common_cells, drop = FALSE]
meta <- meta[match(common_cells, meta$NAME), , drop = FALSE]
coords <- coords[match(common_cells, coords$NAME), , drop = FALSE]
rownames(meta) <- meta$NAME
rownames(coords) <- coords$NAME

sp <- CreateSeuratObject(counts = mat, meta.data = meta)
sp$X <- as.numeric(coords$X)
sp$Y <- as.numeric(coords$Y)
sp$spatial_label <- if ("cluster" %in% colnames(sp@meta.data)) sp@meta.data$cluster else NA_character_

sp <- NormalizeData(sp, verbose = FALSE)

score_mean <- function(seu_obj, genes) {
  genes <- intersect(genes, rownames(seu_obj))
  if (length(genes) == 0) return(rep(NA_real_, ncol(seu_obj)))
  mat <- GetAssayData(seu_obj, slot = "data")
  Matrix::colMeans(mat[genes, , drop = FALSE])
}

safe_scale <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

EA_release_genes <- intersect(c("PLD1", "PLD2", "GDPD1", "GDPD5", "FAAH", "NAPEPLD"), rownames(sp))
EA_recycle_genes <- intersect(c("ETNK1", "ETNK2", "PCYT2", "SELENOI", "CEPT1"), rownames(sp))
FLT1_axis_genes  <- intersect(c("FLT1", "PGF", "KDR", "VEGFA", "ENG", "HIF1A", "EPAS1", "IL1B", "CXCL8", "HMOX1"), rownames(sp))

sp$EA_release <- score_mean(sp, EA_release_genes)
sp$EA_recycle <- score_mean(sp, EA_recycle_genes)
sp$EA_proxy <- safe_scale(sp$EA_release) - safe_scale(sp$EA_recycle)
sp$FLT1_axis <- score_mean(sp, FLT1_axis_genes)

# local-neighborhood smoothing for spatial hotspot visualization
knn_smooth <- function(x, X, Y, k = 15) {
  xy <- cbind(X, Y)
  d <- as.matrix(dist(xy))
  diag(d) <- Inf
  idx <- apply(d, 1, function(z) order(z)[1:k])
  if (is.vector(idx)) idx <- matrix(idx, nrow = 1)
  out <- sapply(seq_along(x), function(i) mean(x[idx[, i]], na.rm = TRUE))
  out
}

sp$EA_proxy_local <- knn_smooth(sp$EA_proxy, sp$X, sp$Y, k = 15)
sp$FLT1_axis_local <- knn_smooth(sp$FLT1_axis, sp$X, sp$Y, k = 15)

plot_spatial <- function(df, val, ttl) {
  ggplot(df, aes(x = X, y = Y, color = .data[[val]])) +
    geom_point(size = 0.6) +
    coord_equal() +
    theme_void() +
    labs(title = ttl, color = val)
}

plot_df <- sp@meta.data %>% mutate(cell = rownames(sp@meta.data))

p1 <- plot_spatial(plot_df, "EA_proxy", "EA proxy")
p2 <- plot_spatial(plot_df, "EA_proxy_local", "EA proxy local hotspot")
p3 <- plot_spatial(plot_df, "FLT1_axis", "FLT1 axis")
p4 <- plot_spatial(plot_df, "FLT1_axis_local", "FLT1 axis local hotspot")

save_plot(file.path(fig_dir, "Fig36_Spatial_EA_FLT1_hotspots.pdf"), (p1 | p2) / (p3 | p4), w = 14, h = 12)

# summary by spatial label
sum_tbl <- plot_df %>%
  group_by(spatial_label) %>%
  summarise(
    n_cells = dplyr::n(),
    EA_proxy = mean(EA_proxy, na.rm = TRUE),
    EA_proxy_local = mean(EA_proxy_local, na.rm = TRUE),
    FLT1_axis = mean(FLT1_axis, na.rm = TRUE),
    FLT1_axis_local = mean(FLT1_axis_local, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(sum_tbl, file.path(tbl_dir, "spatial_EA_FLT1_label_summary.csv"), row.names = FALSE)
qs::qsave(sp, file.path(obj_dir, "architecture_spatial_with_EA_FLT1.qs"))

write_legend(
  "Fig36", "Spatial hotspots of EA proxy and FLT1 axis",
  hypothesis = "Spatially localized EA-rich and vascular-stress niches should exist in first-trimester placental architecture and help explain where Fn-like processes could be most damaging.",
  methods = "We built a Seurat object from the architecture atlas expression and coordinate files, computed EA and FLT1-axis proxy scores, and generated both per-cell and local-neighborhood hotspot maps.",
  readout = "Spatially coherent hotspots near trophoblast/endothelial/fibroblast interfaces would strengthen the niche-based hypothesis.",
  interpretation_template = "- EA hotspots alone suggest metabolic permissiveness.\n- FLT1-axis hotspots alone suggest vascular vulnerability.\n- Overlapping or adjacent hotspots are especially interesting for the toxic-switch model.",
  outfile = file.path(CFG$dirs$legends, "Fig36_legend.md")
)

log_msg("14_spatial_projection_ea_flt1 done.", log_file = log_file)
