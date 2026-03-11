# HooEtAlEnhanced: review + practical run guide

This note summarizes what is usable from this folder, what is mismatched, and how to run the parts that can be trusted in this repository.

## 1) Files are present now

Confirmed present in repo root:
- `1-s2.0-S2405471224001170-main.pdf`
- `Kurkcu_MegL_Toxic_Switch_Hypothesis-Notes.docx`

Confirmed present inside this folder:
- `HooEtAlEnhanced/1-s2.0-S2405471224001170-main.pdf`
- `HooEtAlEnhanced/AdditionalFigures/*.pdf` and `*.png`
- `HooEtAlEnhanced/R-scripts/*.R`

## 2) What looks useful to keep

Useful as **supporting material**:
- `AdditionalFigures/` outputs (good for exploratory/appendix figure ideas).
- `R-scripts/09A_fn_vulnerability_scoring.R`, `09B_organoid_invivo_comparison.R`, `09C_infection_response_analysis.R` as idea templates for added analyses.
- `explant_vs_invivo_celltype_comparison.csv` and `explant_vulnerability_summary.csv` as reference summaries.

## 3) What is mismatched / not directly runnable here

The following scripts are wired to a different project layout and will fail in this repo without path refactoring:
- `R-scripts/RUN_PIPELINE_ENHANCED.R`
- `R-scripts/09A_fn_vulnerability_scoring.R`
- `R-scripts/09B_organoid_invivo_comparison.R`
- `R-scripts/09C_infection_response_analysis.R`

Why:
- they source non-existent files such as `config/config.R`, `scripts/R/utils.R`, and many `scripts/0X_*` paths that are not in this repository.
- several scripts use absolute paths like `/workspace/metadata.csv` and `/workspace/analysis/...`.

## 4) Recommended approach (best practice for your thesis)

1. Keep the canonical repo pipeline as the primary path (`scripts/01` through `scripts/09`).
2. Treat `HooEtAlEnhanced` as a **method/figure idea bank**.
3. Port only specific analysis blocks into the canonical scripts after:
   - replacing absolute paths with `file.path()` based on project root,
   - using your existing `data/seu.qs` object and current metadata harmonization,
   - preserving donor-aware pseudobulk DE on non-integrated counts as the inferential backbone.
4. Use Harmony/RPCA integration only as sensitivity/QC for comparability plots, not as the sole DE basis.

## 5) Minimal run order to stay consistent

From repo root:

```bash
Rscript scripts/01_load_make_seurat.R
Rscript scripts/02_qc_overview.R
Rscript scripts/03_core_umaps_and_composition.R
Rscript scripts/04_gene_sets_scores_plots.R
Rscript scripts/05_susceptibility_severity_models.R
Rscript scripts/06_nk_cytotoxic_module.R
Rscript scripts/07_integration_sensitivity.R
Rscript scripts/08_organoid_vs_placenta_comparison.R
Rscript scripts/09_reproducibility_report.R
```

## 6) If you want to operationalize HooEtAlEnhanced next

Open a follow-up task to do a controlled port:
- add one new script at a time under `scripts/` (not inside `HooEtAlEnhanced`),
- require explicit input checks,
- add outputs under `outputs/` with reproducibility metadata,
- compare key conclusions before/after port to ensure biological consistency.
