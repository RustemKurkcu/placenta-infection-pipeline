# Fig01: QC metrics violin plots

## Hypothesis / idea being tested
We need to verify data quality and identify outliers that could bias downstream analysis (e.g., doublets, dying cells).

## Methods (what we did)
Plot nCount_RNA, nFeature_RNA, and percent.mt across all cells (no filtering yet).

## Readout (what to look for)
Extremely high nCount/nFeature suggests doublets; high percent.mt suggests stressed/dying cells.

## Interpretation template
- If strong tails exist, define filtering thresholds and re-run.
- If QC is reasonable, proceed without aggressive filtering (avoid removing true biology).

