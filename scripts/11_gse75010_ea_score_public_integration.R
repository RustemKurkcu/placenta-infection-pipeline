suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(dplyr)
  library(ggplot2)
  library(qs)
})

# ============================================================
# 11_gse75010_ea_score_public_integration.R
#
# HYPOTHESIS
# ----------
# If the host-side ethanolamine / FLT1 logic is clinically relevant to PE,
# related signatures should be detectable in public placental PE datasets.
# ============================================================

out_dir <- "outputs/public_validation"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gset_list <- getGEO("GSE75010", GSEMatrix = TRUE)
if (length(gset_list) == 0) stop("Could not download GSE75010")
gset <- gset_list[[1]]
expr <- exprs(gset)
pdat <- pData(gset)
fdat <- fData(gset)

if ("gene_assignment" %in% colnames(fdat)) {
  symbols <- sub(" .*", "", fdat$gene_assignment)
} else if ("Gene Symbol" %in% colnames(fdat)) {
  symbols <- fdat[["Gene Symbol"]]
} else if ("Symbol" %in% colnames(fdat)) {
  symbols <- fdat$Symbol
} else {
  stop("Could not identify gene symbol column in feature data")
}

keep <- !is.na(symbols) & symbols != ""
expr <- expr[keep, , drop = FALSE]
symbols <- symbols[keep]
expr_df <- as.data.frame(expr)
expr_df$gene <- symbols
expr_gene <- expr_df %>% group_by(gene) %>% summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")
expr_mat <- as.matrix(expr_gene[, -1])
rownames(expr_mat) <- expr_gene$gene

EA_LIB <- intersect(c("PLD1", "PLD2", "GDPD1", "GDPD5", "FAAH", "NAPEPLD"), rownames(expr_mat))
EA_SINK <- intersect(c("ETNK1", "ETNK2", "PCYT2", "SELENOI", "CEPT1"), rownames(expr_mat))
ANGIO <- intersect(c("FLT1", "PGF", "ENG", "HMOX1", "KDR", "VEGFA", "HIF1A", "EPAS1"), rownames(expr_mat))

zrow <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
score_samples <- function(mat, genes) {
  if (length(genes) == 0) return(rep(NA_real_, ncol(mat)))
  colMeans(mat[genes, , drop = FALSE], na.rm = TRUE)
}

EA_lib_score <- score_samples(expr_mat, EA_LIB)
EA_sink_score <- score_samples(expr_mat, EA_SINK)
EA_balance <- zrow(EA_lib_score) - zrow(EA_sink_score)
ANGIO_score <- score_samples(expr_mat, ANGIO)
FLT1_minus_PGF <- if (all(c("FLT1", "PGF") %in% rownames(expr_mat))) zrow(expr_mat["FLT1", ]) - zrow(expr_mat["PGF", ]) else rep(NA_real_, ncol(expr_mat))

score_df <- data.frame(
  sample = colnames(expr_mat),
  EA_lib_score = EA_lib_score,
  EA_sink_score = EA_sink_score,
  EA_balance = EA_balance,
  ANGIO_score = ANGIO_score,
  FLT1_minus_PGF = FLT1_minus_PGF,
  stringsAsFactors = FALSE
)
score_df <- cbind(score_df, pdat[match(score_df$sample, rownames(pdat)), , drop = FALSE])
write.csv(score_df, file.path(out_dir, "GSE75010_scores.csv"), row.names = FALSE)
qs::qsave(score_df, file.path(out_dir, "GSE75010_scores.qs"), preset = "high")

score_df$group <- NA_character_
for (cc in colnames(score_df)) {
  val <- tolower(as.character(score_df[[cc]]))
  hit_pe <- grepl("preecl|pe", val)
  hit_ctl <- grepl("control|normotens|healthy|normal", val)
  score_df$group[is.na(score_df$group) & hit_pe] <- "PE"
  score_df$group[is.na(score_df$group) & hit_ctl] <- "Control"
}

plot_df <- score_df %>% filter(!is.na(group))
if (nrow(plot_df) > 0) {
  make_box <- function(yvar, ttl, fname) {
    p <- ggplot(plot_df, aes(x = group, y = .data[[yvar]], fill = group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.4) +
      theme_bw() +
      labs(title = ttl, x = NULL, y = yvar)
    ggsave(file.path(out_dir, fname), p, width = 6, height = 5)
  }
  make_box("EA_balance", "GSE75010 EA proxy", "GSE75010_EA_proxy_boxplot.pdf")
  make_box("FLT1_minus_PGF", "GSE75010 FLT1 minus PGF proxy", "GSE75010_FLT1minusPGF_boxplot.pdf")
  make_box("ANGIO_score", "GSE75010 angiogenic stress score", "GSE75010_Angio_score_boxplot.pdf")

  design <- model.matrix(~ 0 + group, data = plot_df)
  colnames(design) <- sub("group", "", colnames(design))
  fit <- lmFit(expr_mat[, plot_df$sample, drop = FALSE], design)
  cont <- makeContrasts(PE_vs_Control = PE - Control, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cont))
  top <- topTable(fit2, number = Inf, sort.by = "P")
  write.csv(top, file.path(out_dir, "GSE75010_PE_vs_Control_limma.csv"))
}

message("Done. Check outputs/public_validation and manually inspect phenotype columns if group assignment needs refinement.")
