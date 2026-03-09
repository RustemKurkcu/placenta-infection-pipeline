# Utilities for clearer Seurat plots

featureplot_white_red <- function(seu, feature, split_by = NULL, title = NULL,
                                  min_cut = "q05", max_cut = "q95",
                                  pt_size = CFG$plotting$pt_size,
                                  raster = TRUE) {
  p <- FeaturePlot(
    seu, features = feature,
    split.by = split_by,
    cols = c("white", "red"),
    min.cutoff = min_cut, max.cutoff = max_cut,
    order = TRUE,
    pt.size = pt_size,
    raster = raster
  )
  if (!is.null(title)) p <- p + ggtitle(title)
  p
}

dimplot_split <- function(seu, group_by, split_by = NULL, title = NULL,
                          pt_size = CFG$plotting$pt_size, raster = TRUE) {
  p <- DimPlot(
    seu, group.by = group_by, split.by = split_by,
    pt.size = pt_size, raster = raster
  )
  if (!is.null(title)) p <- p + ggtitle(title)
  p
}
