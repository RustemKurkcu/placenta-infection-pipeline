# Placenta infection scRNA-seq reanalysis -> Fn organoid planning

## Goal
Reproduce and extend the single-cell explant infection analysis from Hoo et al. (Cell Systems, 2024) as a **reference test-case** for an upcoming *Fusobacterium nucleatum* placental organoid model.

Core questions:
1. **State vs composition:** Do infections change which cell types are present, or mostly change gene programs within the same cell types?
2. **Entry vs response coupling:** Do trophoblast “entry/susceptibility” proxies (adhesion/glycocalyx programs) predict macrophage “severity” (innate/IFN programs)?
3. **Pathogen- and time-specific programs:** What genes and pathways change at **24h vs 48h** for each pathogen, within each key compartment (VCT/EVT/HBC/PAMM)?
4. **Signaling:** Which ligand–receptor axes plausibly connect infected trophoblast states to immune activation?

## What is in this folder
- `01_master_pipeline.R` — end-to-end script that:
  - loads `data/02_processed/seurat_object.qs`
  - writes figures to `results/figures`
  - writes DE tables to `results/tables`
  - writes legends to `results/legends`
  - logs all steps to `results/logs/run.log`

## How to run
From RStudio (recommended):

```r
source("01_master_pipeline.R")
```

Notes:
- The script automatically prefers `stage_perInfection` (matched controls like `UI_Tg_24h`) if present.
- Module scores are computed with **sparse** operations to avoid `R_Calloc` memory errors.
- Pseudobulk DE uses edgeR quasi-likelihood with donor blocking where possible.
