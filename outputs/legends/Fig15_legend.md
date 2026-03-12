# Fig15: QC: NK_like_strict vs others

## Hypothesis / idea being tested
If flagged NK-like cells are doublets, they often show elevated nCount_RNA / nFeature_RNA.

## Methods (what we did)
Compare QC distributions for NK_like_strict vs non-flagged cells.

## Readout (what to look for)
Higher counts/features in flagged group suggests doublets; similar QC supports real cells.

## Interpretation template
- If doublet-like, remove or ignore for biological inference.
- If QC normal and markers coherent, proceed to subset/recluster for activation programs.

