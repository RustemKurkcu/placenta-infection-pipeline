suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(dplyr)
  library(ggplot2)
})

outdir <- "outputs/public_validation"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------
# WHY:
# GSE75010 is a large public placental PE dataset.
# We use it to validate the EA proxy + FLT1 axis in independent data.
# -----------------------------------------
gse <- getGEO("GSE75010", GSEMatrix = TRUE)
eset <- gse[[1]]

expr <- exprs(eset)
pheno <- pData(eset)
feat <- fData(eset)

# inspect phenotype columns once
write.csv(pheno, file.path(outdir, "GSE75010_pheno_raw.csv"))

# try to find a usable gene symbol column
gene_col <- c("Gene Symbol", "GENE_SYMBOL", "Symbol", "gene_assignment")
gene_col <- gene_col[gene_col %in% colnames(feat)][1]
if (is.na(gene_col)) stop("Could not identify gene-symbol column in feature data.")

feat$gene_symbol <- feat[[gene_col]]
keep <- !is.na(feat$gene_symbol) & feat$gene_symbol != ""
expr <- expr[keep, , drop = FALSE]
feat <- feat[keep, , drop = FALSE]

rownames(expr) <- feat$gene_symbol
expr <- avereps(expr)

# crude phenotype parsing; adjust after inspecting GSE75010_pheno_raw.csv
pheno_txt <- apply(pheno, 1, function(x) paste(x, collapse = " | "))
group <- ifelse(grepl("preeclamp", pheno_txt, ignore.case = TRUE), "PE",
                ifelse(grepl("control|normotensive|normal", pheno_txt, ignore.case = TRUE), "Control", NA))
pheno$group_simple <- group

use <- !is.na(pheno$group_simple)
expr <- expr[, use, drop = FALSE]
pheno <- pheno[use, , drop = FALSE]

ea_release <- intersect(c("PLD1","PLD2","GDPD1","GDPD5","FAAH","NAPEPLD"), rownames(expr))
ea_recycle <- intersect(c("ETNK1","ETNK2","PCYT2","SELENOI","CEPT1"), rownames(expr))
flt1_axis  <- intersect(c("FLT1","PGF","KDR","VEGFA","ENG","HIF1A","EPAS1"), rownames(expr))

score_mean <- function(mat, genes) {
  if (length(genes) == 0) return(rep(NA_real_, ncol(mat)))
  colMeans(mat[genes, , drop = FALSE], na.rm = TRUE)
}

scores <- data.frame(
  sample = colnames(expr),
  group = pheno$group_simple,
  EA_release = score_mean(expr, ea_release),
  EA_recycle = score_mean(expr, ea_recycle),
  FLT1_axis  = score_mean(expr, flt1_axis),
  FLT1 = if ("FLT1" %in% rownames(expr)) expr["FLT1", ] else NA_real_,
  PGF  = if ("PGF"  %in% rownames(expr)) expr["PGF", ]  else NA_real_
) %>%
  mutate(
    EA_proxy = EA_release - EA_recycle,
    FLT1_PGF_diff = FLT1 - PGF
  )

write.csv(scores, file.path(outdir, "GSE75010_scores.csv"), row.names = FALSE)

design <- model.matrix(~ 0 + group, data = scores)
colnames(design) <- gsub("^group", "", colnames(design))
fit <- lmFit(expr[, scores$sample], design)
cont <- makeContrasts(PE_vs_Control = PE - Control, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont))
tt <- topTable(fit2, number = Inf, sort.by = "P")
write.csv(tt, file.path(outdir, "GSE75010_PE_vs_Control_limma.csv"))

for (feat in c("EA_proxy", "FLT1_axis", "FLT1", "PGF", "FLT1_PGF_diff")) {
  p <- ggplot(scores, aes(x = group, y = .data[[feat]], fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15) +
    theme_bw() +
    labs(title = paste("GSE75010:", feat))
  ggsave(file.path(outdir, paste0("GSE75010_", feat, ".pdf")), p, width = 6, height = 5)
}