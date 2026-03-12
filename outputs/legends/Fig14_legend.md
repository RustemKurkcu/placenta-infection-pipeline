# Fig14: UMAP: NK_like_strict flagged cells

## Hypothesis / idea being tested
If infection recruits/expands NK-like cells or induces NK-like cytotoxic states, flagged cells should increase by infection/time and form a coherent cluster.

## Methods (what we did)
Compute strict NK-like boolean gate from expression fetched via FetchData; plot flagged cells on UMAP; tabulate by condition.

## Readout (what to look for)
A coherent cluster + condition-associated increase supports a real signal; scattered rare cells suggests contamination/doublets.

## Interpretation template
- If extremely rare (<0.1%) and scattered, treat as likely contamination.
- If enriched in a condition, re-cluster immune cells and validate with additional markers (FCGR3A, TRBC1/2, TRAC, MS4A7, etc.).

