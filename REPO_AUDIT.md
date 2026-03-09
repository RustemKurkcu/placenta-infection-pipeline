# Repository audit (pipeline, stats, and biology)

## Executive assessment

- The repository is **partly updated well** for analysis code structure: there is a modular 6-step pipeline with shared config/helpers and explicit output paths under `outputs/`.
- The repo is **not yet fully reproducible** in its current committed state because key runtime inputs are missing from versioned tree (e.g., `data/*.qs` objects and `gene_sets/*.txt` files referenced by scripts).
- There are also **two competing pipelines** (`scripts/01_master_pipeline.R` vs. modular `scripts/01..06` + `R/` helpers), which risks divergence and inconsistent results.

## What is good

1. Metadata harmonization is explicit in load step (`hpi_num`, `condition` creation), which is needed for downstream stratification.
2. Pseudobulk DE in the master script is conceptually strong: aggregate by replicate, include donor in design, use edgeR QL robust fitting.
3. LFS patterns are configured for `.qs` and `.rds` files.

## Main issues to fix before trusting final conclusions

1. **Missing required inputs in repo**
   - `scripts/04_gene_sets_scores_plots.R` expects:
     - `gene_sets/glyco_genes.txt`
     - `gene_sets/adhesion_genes.txt`
     - `gene_sets/innate_genes.txt`
   - The committed tree currently does not include `gene_sets/`.

2. **Path mismatch between pipelines**
   - `scripts/01_master_pipeline.R` reads `data/02_processed/seurat_object.qs` and writes to `results/`.
   - Modular scripts read/write `data/` + `outputs/`.
   - This can lead to people running different workflows and obtaining incomparable outputs.

3. **Legacy monolithic script has a logic bug**
   - In `scripts/placenta_fn_pipeline.R`, `subset(seu_obj, subset = cell_type == cell_type)` will always be TRUE for non-missing values, so DE is not actually cell-type specific.

4. **Statistical fragility in current sample-level correlation endpoint**
   - In this snapshot, several infection×time strata appear to have only 2–4 samples, which is insufficient for stable correlations.
   - Report these as exploratory, not confirmatory.

## Statistical and biological validity verdict (current snapshot)

- **QC distributions**: broadly plausible for scRNA-seq (median nCount/nFeature and mito proportions are in expected ranges), but this alone does not rule out doublets/ambient RNA.
- **Susceptibility→severity model**: biologically plausible framing, but currently vulnerable to confounding unless analyses are per-pathogen/time and donor-aware.
- **NK-like conclusion**: strict NK-like gate is extremely rare in this dataset snapshot (very low fraction), so strong biological claims about NK expansion are not yet supported.

## Minimal additions to make this publishable/reproducible

1. Add `gene_sets/` text files used by script 04.
2. Add a top-level `README.md` with:
   - canonical entrypoint (choose one pipeline),
   - required input object schema,
   - exact command sequence,
   - expected output file list.
3. Add a machine-readable run manifest (`outputs/logs/session_info.txt` + package versions).
4. Add one summary table with per-condition sample counts (donors and pseudobulk replicates) used in each statistical test.

## Recommendation

- Treat `scripts/placenta_fn_pipeline.R` as legacy and run the modular pipeline only (or archive legacy script).
- Keep the current biological conclusions at the level of **hypothesis-generating** until replicate counts are strengthened and DE/correlation results are cross-validated.
