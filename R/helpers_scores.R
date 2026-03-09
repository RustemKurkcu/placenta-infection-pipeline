present_genes <- function(seu, genes) {
  genes <- unique(genes)
  genes[genes %in% rownames(seu)]
}

# Low-memory signature scoring: mean expression of gene set in normalized data.
# This is NOT the same as Seurat's AddModuleScore (which adds matched controls).
add_signature_mean <- function(seu, genes, name) {
  genes <- present_genes(seu, genes)
  if (length(genes) == 0) stop("No genes found for signature: ", name)

  mat <- tryCatch(
    GetAssayData(seu, layer = "data"),
    error = function(e) GetAssayData(seu, slot = "data")
  )
  score <- Matrix::colMeans(mat[genes, , drop = FALSE])
  seu[[name]] <- score
  seu
}

# More Seurat-like module score, but restrict pool to HVGs to reduce memory.
# Falls back to mean-score if AddModuleScore fails.
add_signature_hvg_pool <- function(seu, genes, name, n_hvg = 2000, ctrl = 50) {
  genes <- present_genes(seu, genes)
  if (length(genes) == 0) stop("No genes found for signature: ", name)

  if (length(VariableFeatures(seu)) < 50) {
    seu <- FindVariableFeatures(seu, nfeatures = n_hvg, verbose = FALSE)
  }
  pool <- VariableFeatures(seu)

  out <- tryCatch({
    seu2 <- AddModuleScore(seu, features = list(genes), name = name, pool = pool, ctrl = ctrl, verbose = FALSE)
    col <- paste0(name, "1")
    names(seu2@meta.data)[names(seu2@meta.data) == col] <- name
    seu2
  }, error = function(e) {
    message("AddModuleScore failed for ", name, ": ", conditionMessage(e))
    message("Falling back to mean signature scoring (no matched controls).")
    add_signature_mean(seu, genes, name)
  })

  out
}

zscore_within <- function(seu, score_col, group_col, new_col = NULL) {
  if (is.null(new_col)) new_col <- paste0(score_col, "_z_", group_col)
  md <- seu@meta.data
  md[[new_col]] <- ave(md[[score_col]], md[[group_col]], FUN = function(x) as.numeric(scale(x)))
  seu@meta.data <- md
  seu
}
