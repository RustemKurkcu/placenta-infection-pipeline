# Placenta infection pipeline

Canonical workflow (modular):

1. `scripts/01_load_make_seurat.R`
2. `scripts/02_qc_overview.R`
3. `scripts/03_core_umaps_and_composition.R`
4. `scripts/04_gene_sets_scores_plots.R`
5. `scripts/05_susceptibility_severity_models.R`
6. `scripts/06_nk_cytotoxic_module.R`
7. `scripts/07_integration_sensitivity.R`
8. `scripts/08_organoid_vs_placenta_comparison.R`
9. `scripts/09_reproducibility_report.R`

## Required inputs

- `data/seu.qs` (or `data/seu.rds`, configurable in `R/config.R`)
- gene sets:
  - `gene_sets/glyco_genes.txt`
  - `gene_sets/adhesion_genes.txt`
  - `gene_sets/innate_genes.txt`

## Notes

- Use non-integrated RNA counts for pseudobulk DE inference.
- Use integration (RPCA/Harmony) as a sensitivity analysis for comparability; see `INTEGRATION_COMPARABILITY_GUIDE.md`.
- `scripts/placenta_fn_pipeline.R` is retained as a legacy compatibility script.


## Advanced guidance

- Integration/comparability: `INTEGRATION_COMPARABILITY_GUIDE.md`
- Enhancement roadmap: `ANALYSIS_ENHANCEMENT_PLAN.md`
