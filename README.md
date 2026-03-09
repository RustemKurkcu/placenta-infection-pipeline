# Placental infection scRNA-seq pipeline (Seurat)

This project is a **reproducible, organized** workflow to:
1) load a Seurat object (or build one from h5ad / CellxGene downloads),
2) generate core figures (UMAPs, composition),
3) compute **proxy signatures** (glycocalyx/glyco, adhesion, innate severity),
4) test the **susceptibility vs severity** idea at cell and pseudo-bulk (sample) levels,
5) audit the **"NK-like / cytotoxic-like"** signal and decide whether it is real or artifact.

## Quick start

1. Open RStudio in this project folder (or set working directory here).
2. Install packages:
   - run `00_setup/00_install_packages.R`
3. Run the pipeline:
   - `source("scripts/01_load_make_seurat.R")`
   - `source("scripts/02_qc_overview.R")`
   - `source("scripts/03_core_umaps_and_composition.R")`
   - `source("scripts/04_gene_sets_scores_plots.R")`
   - `source("scripts/05_susceptibility_severity_models.R")`
   - `source("scripts/06_nk_cytotoxic_module.R")`

Outputs are written to `outputs/`:
- `outputs/figures/`  (PDF/PNG figures)
- `outputs/legends/`  (one Markdown file per figure: hypothesis, methods, readout, interpretation template)
- `outputs/tables/`   (CSVs)
- `outputs/objects/`  (saved Seurat objects)

## Where to put your data

- Put your Seurat object in `data/` (recommended name: `seu.rds` or `seu.qs`)
- If starting from h5ad, edit `scripts/01_load_make_seurat.R` and set `CFG$data$h5ad_path`.

## Core hypothesis

**Susceptibility hypothesis:** cell states with higher *entry/attachment surface programs*
(e.g., epithelial adhesion, glycocalyx / glycan / proteoglycan features) are **more permissive**
to placental pathogen interaction.

**Severity hypothesis:** infections that invade/replicate more effectively trigger stronger
**innate inflammatory programs** (in immune and/or trophoblast compartments) and may shift
cell-type composition.

The scripts operationalize these ideas with:
- **proxy gene sets** (glyco/adhesion/innate),
- per-cell signature overlays and per-sample pseudo-bulk summaries,
- comparisons across **infection (UI/Lm/Pf/Tg)** and **time (24h/48h)**.
