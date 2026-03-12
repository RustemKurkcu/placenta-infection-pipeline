# Fig02: QC by condition

## Hypothesis / idea being tested
Some infections/time points may alter RNA content or stress signatures; we want to rule out technical artifacts.

## Methods (what we did)
Violin plots of QC metrics grouped by condition (infection × time).

## Readout (what to look for)
Systematic shifts in QC can confound differential expression and signature scoring.

## Interpretation template
- If one condition has much higher counts/mito, consider per-condition normalization or sensitivity analyses.
- If stable, comparisons are more trustworthy.

