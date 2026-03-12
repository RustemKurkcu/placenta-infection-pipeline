# Script review and next steps (11–15 focus)

## Scope checked in this repository snapshot

- `scripts/01` to `scripts/09` exist.
- `scripts/13_fib2_reference_mapping.R` exists.
- `scripts/14_make_readable_embedding_figures.R` and `scripts/15_figure_readability_audit.R` were added to make 14/15 operational.
- `scripts/11_*` and `scripts/12_*` are not present in this checkout.

## About `concersatonswithchatforupdatingapproach.txt`

That file is not present in this repository snapshot, so ideas in that text file could not be reviewed directly here.

## Evaluation of current approach

### Strong parts

1. **Core pipeline modularity (01–09)** is clear and reusable.
2. **Integration sensitivity (07)** now has Harmony API compatibility plus speed/memory safeguards.
3. **Architecture transfer (13)** uses conservative refinement rules and keeps author labels by default.
4. **Runner tooling (`local_runner`)** provides reproducible execution and output collection.

### Risks / caveats

1. **11 and 12 missing**: if your plan expects numbered stages 11/12, this repo currently has a numbering gap.
2. **Reference dependency for 13**: script 13 requires a Slide-tags reference object in expected paths.
3. **Visualization readability**: earlier figure sets can be hard to interpret when panel density is high or dimensions are too small.

## 14 and 15 status

### 14 (`scripts/14_make_readable_embedding_figures.R`)

Purpose:
- Rebuild embedding plots with larger canvases/themes for human readability.
- Uses best available object (`seu_with_architecture_transfer` → `seu_with_integration_checks` → `seu_with_scores` → `seu_clean`).
- Writes both PDF and PNG outputs.

Expected outputs:
- `outputs/figures/FigReadable_*`
- `outputs/legends/FigReadable_legend.md`

### 15 (`scripts/15_figure_readability_audit.R`)

Purpose:
- Audit figure outputs and flag likely low-detail files + missing legends.

Expected outputs:
- `outputs/tables/figure_readability_audit.csv`
- `outputs/tables/figure_readability_flagged.csv`
- `outputs/tables/figure_readability_audit.md`

## Recommended run order now

1. Canonical: `01..09`
2. Optional transfer refinement: `13`
3. Readability regeneration: `14`
4. Readability audit: `15`

## Results interpretation so far (high-level)

- Current outputs support a **workable exploratory workflow** with cautious integration and conservative label transfer.
- Claims should still be presented as **confidence-weighted** where sample strata are small.
- Use script 15 flags to identify specific figures to reformat before thesis/manuscript finalization.

## Concrete next improvements

1. Add explicit `scripts/11_*` and `scripts/12_*` (or renumber docs) to remove numbering ambiguity.
2. Add a run summary table (n cells, n donors, % refined labels, % cross-lineage overrides).
3. Add fixed theme/palette standards to all figure-producing scripts for consistency.
