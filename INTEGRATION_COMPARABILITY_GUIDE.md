# Integration and cross-sample comparability guide

## What to use here: integration vs. no integration

For this placenta infection design (multiple donors, infection states, timepoints), a good default is:

1. Keep **RNA assay (non-integrated)** for DE/pseudobulk statistics.
2. Use an integrated latent space only for visualization / clustering robustness checks.

This avoids over-correcting true infection biology.

## Recommended workflow

1. **Primary analysis (already aligned with your scripts)**
   - Use donor-aware pseudobulk edgeR on raw counts (current approach in `01_master_pipeline.R`).
   - Compare matched controls (`stage_perInfection` style) whenever available.

2. **Integration sensitivity analysis (add as QC, not replacement)**
   - Build one integrated embedding with Seurat RPCA or Harmony using donor as batch.
   - Re-check whether major cell-type structure is preserved and whether infection effects remain visible within cell types.

3. **Decide if stronger integration is needed**
   - If donor dominates PCA/UMAP and masks lineage structure: try RPCA/Harmony.
   - If datasets are highly asymmetric or multi-study: consider scVI (outside current repo scope).

## Method choice quick matrix

- **Seurat RPCA integration**: good first choice when datasets are similar and you want conservative correction.
- **Harmony**: convenient when you already have PCA and want fast donor-effect correction.
- **CCA integration**: older Seurat path; can work, but RPCA is usually preferred for speed/stability.
- **scVI**: strong for complex batch effects but adds Python dependency and modeling complexity.

## Acceptance checks before cross-sample comparisons

Run these checks before final claims:

1. Per-cell-type donor mixing in embedding improves **without collapsing infection/time separation**.
2. Marker genes for known cell types remain specific after integration.
3. Pseudobulk DE direction/significance is qualitatively stable (integrated-vs-nonintegrated clustering assignments).
4. Per-contrast replicate count table is reported (donor count and pseudobulk replicates).
5. Any strata with n <= 2 are labeled exploratory only.

## Bottom line for your current repo

- Your current donor-aware pseudobulk strategy is statistically appropriate as the main inferential engine.
- Add one integration sensitivity pass (RPCA or Harmony) to demonstrate conclusions are not driven by donor structure.
