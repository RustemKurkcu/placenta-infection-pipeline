# Fig35 — Architecture informed label refinement

**Panel A — Author labels.** UMAP projection of the infection dataset colored by the original author-supplied cell types (`celltype_author`).

**Panel B — Refined labels after architecture transfer.** The same UMAP colored by `celltype_refined` after conservative architecture-based label transfer from the Slide-tags mapped multiome reference. Refinement rules applied: author labels are retained by default; within-lineage refinements occur at `score >= 0.85`; cross-lineage overrides allowed only at `score >= 0.95` and `margin >= 0.10`.

**Panel C — Transfer confidence.** Continuous color scale displays maximum transfer confidence (`arch_prediction_score_max`).

**Interpretation notes:**
- Concentration of FIB2 predictions within stromal author labels supports the fibroblast reservoir hypothesis.
- Consult `architecture_transfer_top3_predictions.csv` for per-cell audit.
