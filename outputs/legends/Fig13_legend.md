# Fig13: Marker panel for NK/T/myeloid/trophoblast identity

## Hypothesis / idea being tested
If true NK/T cells are present, they should show coherent immune identity (PTPRC+) and lineage markers (TRAC/CD3D for T; NKG7/KLRD1/GNLY/PRF1 for NK) without trophoblast contamination (KRT8/18).

## Methods (what we did)
Overlay immune and cytotoxic markers on the global UMAP; check co-localization and exclusivity.

## Readout (what to look for)
A real lymphocyte cluster would be PTPRC+ and TRAC/CD3D+ (T) or PTPRC+ and NKG7/KLRD1/GNLY+ (NK), and largely KRT-.

## Interpretation template
- If markers are sparse and not co-expressed, likely ambient RNA or doublets.
- If a coherent cluster exists, subset and re-cluster that compartment for deeper profiling.

