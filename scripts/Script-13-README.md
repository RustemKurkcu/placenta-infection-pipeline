# Script 13 README — Architecture-informed label transfer + postprocessing

## Overview

T core Script **`13_fib2_reference_mapping.R`** and a companion post-processing script (e.g. `13_postprocess_figures_and_tables.R`) that produce refined cell-type annotations for the infection dataset using the Slide-tags placenta atlas as a high-resolution reference, plus a set of audit tables and figures to validate the refinements.

**Primary scientific goal / hypothesis**

- Test whether a **FIB2-like stromal subpopulation** (defined in the first-trimester placenta atlas) exists within the infection samples’ stromal/perivascular compartments.
- Produce **conservative**, auditable refinements of author/paper cell-type labels using architecture-based label transfer — only apply refinements when transfer confidence is high and lineage-consistent.
- Provide outputs (per-cell top-3 predictions and scores, confusion matrices, per-author summaries, marker tables and figures) that fully document how refinements were made.

---

## Scripts

### `13_fib2_reference_mapping.R`
**Purpose:** perform label transfer (Seurat anchors) from Slide-tags multiome reference to the infection dataset, record top-3 predictions and scores per cell, and apply conservative refinement rules to produce `celltype_refined` and `fib_subtype`.

**Main steps:**
1. Load repo config (`R/config.R`) and helper functions (`R/helpers_io.R`); load Seurat, dplyr, ggplot2, qs, patchwork, (SingleR optional).
2. Choose a query object in `outputs/objects/` (priority list: `seu_with_ea_flt1_proxies.qs`, `seu_with_integration_checks.qs`, `parameter_sweep/seu_HVG4000_PC40.qs`, `seu_HARM.qs`). Merge proxy metadata when appropriate.
3. Load the Slide-tags reference (`data/02_processed/slidetags_mapped_to_multiome.rds` or equivalent).
4. Intersect features between query and reference (require ≥ 500 shared genes) and optionally subset query to relevant lineages (stromal, EVT, VCT, immune) for memory efficiency.
5. Compute PCA in reference and query (only if missing) via `prep_for_transfer`.
6. Build transfer anchors: `FindTransferAnchors(reference=seu_ref_sub, query=seu_inf_sub, dims=1:30, normalization.method="LogNormalize")`.
7. Predict labels via `TransferData(..., k.weight=50, weight.reduction="pcaproject")`. Extract per-label scores and `prediction.score.max`.
8. Compute **top-3 predicted labels** and **margin** (`score1 - score2`) per cell.
9. Map predictions back to the full `seu_inf` object, and apply **conservative refinement rules**:
   - Always retain `celltype_author` by default.
   - Fill missing author labels if `score1 >= PRED_SCORE_HIGH` (default 0.70).
   - Within-lineage refinement if `author_lineage == pred_lineage` and `score1 >= PRED_SCORE_REFINE` (default 0.85).
   - Cross-lineage override only if `score1 >= CROSS_LINEAGE_SCORE` (default 0.95) **and** `(score1 - score2) >= CROSS_LINEAGE_MARGIN` (default 0.10).
10. Add `fib_subtype` for FIB1/FIB2 labels and optional immune subclassification (SingleR + module scores) when an immune reference is available.
11. Write audit CSVs, figures, and a canonical Seurat QS (`seu_with_architecture_transfer.qs`).

**Key inputs**
- Query (infection) object (one of):  
  `outputs/objects/seu_with_ea_flt1_proxies.qs` (preferred),  
  `outputs/objects/seu_HVG4000_PC40.qs`, or similar.
- Reference object:  
  `data/02_processed/slidetags_mapped_to_multiome.rds` (expected to contain `cell_type_cluster` and `predicted.id`).
- Optional immune reference: `data/processed/immune_reference.qs`.

**Key outputs**
Stored under `outputs/tables/`, `outputs/figures/`, and `outputs/objects/`:
- `architecture_transfer_summary.csv` — counts by `(celltype_author, celltype_refined, refinement_action)`.
- `architecture_transfer_prediction_scores.csv` — per-refined-label mean/median scores.
- `architecture_transfer_top3_predictions.csv` — per-cell `pred1/score1`, `pred2/score2`, `pred3/score3`, `margin`.
- `architecture_transfer_threshold_sensitivity.csv` — sensitivity of refinements to thresholds.
- `refinement_action_by_author.csv` — breakdown of refinement actions per author label.
- `Fig35_Architecture_transfer_comparison.pdf` — side-by-side UMAPs: author | refined | transfer confidence.
- `Fig35b_Architecture_transfer_confidence.pdf` — transfer-confidence UMAP.
- `outputs/objects/seu_with_architecture_transfer.qs` — full Seurat object with new metadata.

**Parameters & recommended defaults**
- `DIMS = 1:30` (PCA dims for anchors; lower for low memory)
- `K.WEIGHT = 50` (k weight for TransferData)
- `PRED_SCORE_HIGH = 0.70` (fill missing author labels)
- `PRED_SCORE_REFINE = 0.85` (within-lineage refinement)
- `CROSS_LINEAGE_SCORE = 0.95`, `CROSS_LINEAGE_MARGIN = 0.10` (strict cross-lineage override)
- `RUN_IMMUNE_SUBCLASS = TRUE` (optional SingleR immune refinement)

**Rationale**
- The Slide-tags multiome reference has biologically meaningful fibroblast (FIB1/FIB2) and immune labels. However, to avoid mislabeling, we adopt a conservative, lineage-aware refinement policy: only accept predictions when they are high-confidence and lineage-consistent.

---

### `13_postprocess_figures_and_tables.R` (companion post-processing)

**Purpose:** produce extra audit tables and publication-ready figures that quantitatively and visually validate the transfer and refinements.

**Main outputs produced:**
- `confusion_author_to_refined_counts.csv` — raw counts author → refined.
- `confusion_author_to_refined_pct_by_author.csv` — normalized per-author (row-wise).
- `Fig36_confusion_heatmap.png` — heatmap of the confusion matrix (proportion by author).
- `per_author_summary.csv` — per-author `n_cells`, `n_refined`, `prop_refined`, score summary, cross-lineage override counts.
- `Fig37_prediction_score_by_author.pdf` — violin/boxplot of transfer confidence across author types.
- `Fig38_refinement_action_by_author.pdf` — stacked barplots of refinement actions per author label.
- `architecture_transfer_top3_predictions.csv` (regenerated if missing).
- `architecture_transfer_threshold_sensitivity.csv` (regenerated).
- `FIB2_vs_FIB1_markers.csv` or `markers_by_refined_label.csv` — differential markers for refined labels.
- `Fig39_marker_DotPlot_by_refined.pdf` — DotPlot of canonical markers across refined labels.
- `Fig40_fibro_featureplots.pdf` — UMAP featureplots for fibro markers (PDGFRA, DCN, CXCL14, SPON1) for stromal subset.
- `Fig35_legend.md` and `scripts/13_mapping/13_methods.md` — enhanced figure legend & methods description.

**Why these are useful**
- The confusion matrix and per-author summaries quantify how often author labels and refined labels agree and explain any systematic disagreements.  
- Violin/boxplots of confidence per author highlight which author labels are predicted robustly.  
- Marker tables (FindMarkers) and DotPlots/FeaturePlots provide biological validation of refined labels (e.g., PDGFRA/CXCL14 enriched in FIB2).

---

## Figure legends (for manuscript / README)

### Fig35 — Architecture-informed label refinement (three panels)
**Panel A (left):** UMAP of infection dataset colored by author labels (`celltype_author`).  
**Panel B (middle):** Same UMAP colored by `celltype_refined` after architecture transfer; only conservative, high-confidence refinements are applied.  
**Panel C (right):** UMAP colored by maximum transfer confidence (`arch_prediction_score_max`), showing spatial distribution of high/low-confidence predictions.  
**Interpretation:** Concordance between A and B, and concentration of FIB2 predictions in stromal author labels, supports the FIB2 reservoir hypothesis. Consult `architecture_transfer_top3_predictions.csv` for per-cell audit.

---

## How to run (example commands)

From the repository root:

```bash
# run the transfer + refinement
Rscript scripts/13_mapping/13_fib2_reference_mapping.R

# run the postprocessing (tables + extra figures)
Rscript scripts/13_mapping/13_postprocess_figures_and_tables.R