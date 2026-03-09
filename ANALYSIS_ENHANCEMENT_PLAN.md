# Analysis enhancement plan (paper/thesis aligned)

## Important note on requested source documents

I searched the repository snapshot for:
- `1-s2.0-S2405471224001170-main.pdf`
- `Kurkcu_MegL_Toxic_Switch_Hypothesis-Notes.docx`

They are not present in this checkout, so this plan is based on your current code/data files only.

## What was improved now

1. Added integration sensitivity script (`scripts/07_integration_sensitivity.R`) to evaluate donor correction with Harmony while preserving infection biology.
2. Added organoid-vs-placenta mapping script (`scripts/08_organoid_vs_placenta_comparison.R`) using anchor transfer and joint UMAP.
3. Added reproducibility/reporting script (`scripts/09_reproducibility_report.R`) to emit donor/condition replicate tables and `session_info.txt`.

## Nature-grade figure package to target

1. Study design + cohort balance figure
   - Donors/conditions/timepoints and cell counts.
2. Reference atlas overview
   - UMAP + marker validation panels.
3. Infection-state shifts within major lineages
   - Per-lineage pseudobulk volcano + top pathways.
4. Susceptibility->severity synthesis
   - Effect-size forest plot across pathogens/timepoints.
5. Organoid validation against placenta
   - Joint embedding, label-transfer confidence, marker concordance heatmap.

## Statistical upgrades recommended

1. Report effect sizes + confidence intervals (not just p-values).
2. Use per-contrast minimum replicate thresholds and mark underpowered contrasts as exploratory.
3. Add sensitivity analysis:
   - non-integrated baseline
   - donor-corrected latent space (Harmony/RPCA)
   - consistency check of direction/significance.
4. Add pathway-level interpretation for DE genes (camera/fgsea style).

## Organoid vs placenta comparison checklist

1. Harmonize gene symbols and QC thresholds across datasets.
2. Map organoid to placenta reference with `FindTransferAnchors` / `TransferData`.
3. Quantify mapping confidence and ambiguous states.
4. Compare key signatures (glyco, adhesion, innate) in matched predicted cell types.
5. Validate whether model captures expected trophoblast compartments and inflammatory response range.

## Inputs still needed from you

1. Add the PDF and thesis DOCX into repo (or paste key methods/claims sections).
2. Provide organoid Seurat object at `data/organoid_seu.qs`.
3. Confirm final figure priority (main text vs supplement shortlist).
