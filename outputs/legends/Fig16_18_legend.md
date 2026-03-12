# Fig16-18: Integration sensitivity with Harmony

## Hypothesis / idea being tested
Cross-donor comparability can improve after batch correction without erasing infection biology.

## Methods (what we did)
Run Harmony on PCA using donor as batch variable; compare cell-type, donor, and infection overlays.

## Readout (what to look for)
Desired pattern: reduced donor-driven segregation while preserving cell-type structure and condition-associated shifts.

## Interpretation template
- If donor separation dominates after correction, revisit QC and batch variables.
- If infection signal disappears globally, correction may be too aggressive.

