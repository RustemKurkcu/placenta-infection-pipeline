source("R/config.R")
source("R/helpers_io.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(qs)
})

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("08_organoid_vs_placenta_comparison starting.", log_file = log_file)

# Expected inputs:
# - Placenta reference: outputs/objects/seu_with_scores.qs
# - Organoid query: data/organoid_seu.qs (same gene naming convention)
placenta_path <- file.path(CFG$dirs$objects, "seu_with_scores.qs")
organoid_path <- file.path(CFG$dirs$data, "organoid_seu.qs")

if (!file.exists(placenta_path)) stop("Missing placenta reference object: ", placenta_path)
if (!file.exists(organoid_path)) stop("Missing organoid object: ", organoid_path)

ref <- qs::qread(placenta_path)
qry <- qs::qread(organoid_path)
DefaultAssay(ref) <- "RNA"
DefaultAssay(qry) <- "RNA"

ref$dataset <- "placenta"
qry$dataset <- "organoid"

# Keep shared genes for mapping consistency
shared <- intersect(rownames(ref), rownames(qry))
ref <- subset(ref, features = shared)
qry <- subset(qry, features = shared)

# Normalize and PCA on reference
ref <- NormalizeData(ref, verbose = FALSE)
ref <- FindVariableFeatures(ref, nfeatures = 3000, verbose = FALSE)
ref <- ScaleData(ref, features = VariableFeatures(ref), verbose = FALSE)
ref <- RunPCA(ref, features = VariableFeatures(ref), npcs = 50, verbose = FALSE)
ref <- RunUMAP(ref, dims = 1:30, reduction = "pca", reduction.name = "umap_ref", verbose = FALSE)

# Anchor mapping (reference -> query)
anchors <- FindTransferAnchors(
  reference = ref,
  query = qry,
  dims = 1:30,
  reference.reduction = "pca",
  normalization.method = "LogNormalize"
)

pred <- TransferData(
  anchorset = anchors,
  refdata = ref[[CFG$cols$cell_type]][,1],
  dims = 1:30
)
qry <- AddMetaData(qry, metadata = pred)

# Build merged object for side-by-side visual checks
merged <- merge(ref, qry)
merged <- NormalizeData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, nfeatures = 3000, verbose = FALSE)
merged <- ScaleData(merged, features = VariableFeatures(merged), verbose = FALSE)
merged <- RunPCA(merged, features = VariableFeatures(merged), npcs = 50, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1:30, reduction = "pca", reduction.name = "umap_joint", verbose = FALSE)

p1 <- DimPlot(merged, reduction = "umap_joint", group.by = "dataset", raster = TRUE) +
  ggtitle("Joint UMAP: placenta vs organoid")
p2 <- DimPlot(merged, reduction = "umap_joint", group.by = CFG$cols$cell_type, raster = TRUE, label = TRUE, repel = TRUE) +
  ggtitle("Joint UMAP: reference cell type labels")

save_plot(file.path(CFG$dirs$figures, "Fig19_Joint_UMAP_dataset.pdf"), p1, w = 12, h = 8)
save_plot(file.path(CFG$dirs$figures, "Fig20_Joint_UMAP_celltype.pdf"), p2, w = 12, h = 8)

# Confidence table for transferred labels
if ("prediction.score.max" %in% colnames(qry@meta.data)) {
  tab <- qry@meta.data %>%
    count(predicted.id, wt = prediction.score.max, name = "sum_prediction_score") %>%
    arrange(desc(sum_prediction_score))
  ensure_dir(CFG$dirs$tables)
  write.csv(tab, file.path(CFG$dirs$tables, "organoid_label_transfer_scores.csv"), row.names = FALSE)
}

qs::qsave(qry, file.path(CFG$dirs$objects, "organoid_mapped_to_placenta.qs"))
qs::qsave(merged, file.path(CFG$dirs$objects, "placenta_organoid_joint.qs"))

write_legend(
  "Fig19-20", "Organoid-to-placenta mapping and joint embedding",
  hypothesis = "If organoid model recapitulates placental states, organoid cells should map to plausible placenta cell types with reasonable confidence.",
  methods = "Reference mapping with FindTransferAnchors/TransferData and joint UMAP visualization on shared genes.",
  readout = "Look for organoid overlap with expected trophoblast/immune/stromal compartments and inspect prediction confidence.",
  interpretation_template = "- High-confidence mapping to expected trophoblast states supports model validity.\n- Systematic shifts indicate model-specific states or missing microenvironmental cues.",
  outfile = file.path(CFG$dirs$legends, "Fig19_20_legend.md")
)

log_msg("08_organoid_vs_placenta_comparison done.", log_file = log_file)
