# Fig03: UMAP by cell type

## Hypothesis / idea being tested
Placental compartments (trophoblast, macrophage, stromal, endothelial) segregate into distinct transcriptomic states.

## Methods (what we did)
Plot existing UMAP embedding and color by the provided cell_type annotation.

## Readout (what to look for)
Well-separated, interpretable clusters suggest good annotation and biology; mixing may indicate batch effects or ambiguous states.

## Interpretation template
- If cell types are well separated, proceed with within-cell-type comparisons.
- If not, consider re-clustering or checking marker genes.

