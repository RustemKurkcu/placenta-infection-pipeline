# Repository audit (pipeline, stats, and biology)

## Executive assessment

- The repository is **partly updated well** for analysis code structure: there is a modular 6-step pipeline with shared config/helpers and explicit output paths under `outputs/`.
- The repo is **not yet fully reproducible** in its current committed state because large runtime inputs (`data/*.qs`) are not versioned in Git and must be provided externally.
- There are also **two competing pipelines** (`scripts/01_master_pipeline.R` vs. modular `scripts/01..06` + `R/` helpers), which risks divergence and inconsistent results.

## What is good

1. Metadata harmonization is explicit in load step (`hpi_num`, `condition` creation), which is needed for downstream stratification.
2. Pseudobulk DE in the master script is conceptually strong: aggregate by replicate, include donor in design, use edgeR QL robust fitting.
3. LFS patterns are configured for `.qs` and `.rds` files.

## Main issues to fix before trusting final conclusions

1. **Statistical fragility in current sample-level correlation endpoint**
   - In this snapshot, several infection×time strata appear to have only 2–4 samples, which is insufficient for stable correlations.
   - Report these as exploratory, not confirmatory.

## Statistical and biological validity verdict (current snapshot)

- **QC distributions**: broadly plausible for scRNA-seq (median nCount/nFeature and mito proportions are in expected ranges), but this alone does not rule out doublets/ambient RNA.
- **Susceptibility→severity model**: biologically plausible framing, but currently vulnerable to confounding unless analyses are per-pathogen/time and donor-aware.
- **NK-like conclusion**: strict NK-like gate is extremely rare in this dataset snapshot (very low fraction), so strong biological claims about NK expansion are not yet supported.

## Minimal additions to make this publishable/reproducible

1. Add a machine-readable run manifest (`outputs/logs/session_info.txt` + package versions).
2. Add one summary table with per-condition sample counts (donors and pseudobulk replicates) used in each statistical test.
3. Add one integration sensitivity notebook/report (RPCA or Harmony) that demonstrates donor correction does not erase infection biology.

## Recommendation

- Treat `scripts/placenta_fn_pipeline.R` as legacy and run the modular pipeline only (or archive legacy script).
- Keep the current biological conclusions at the level of **hypothesis-generating** until replicate counts are strengthened and DE/correlation results are cross-validated.
