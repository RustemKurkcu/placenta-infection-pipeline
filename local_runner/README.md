# Local Runner Bundle (manifest-aware)

This folder gives you a clean, reproducible way to run the pipeline on your own workstation and collect outputs in one logical place.

## What this bundle does

1. **Prepare inputs** from the files you already have (large `.qs`/`.rds` objects).
2. **Run canonical scripts in order** (`scripts/01..09`) with per-step logs.
3. **Collect figures/tables/logs** from both `outputs/` and `results/` into one timestamped run folder.

## Recommended run order

From repo root:

```bash
bash local_runner/01_prepare_inputs.sh
bash local_runner/02_run_pipeline.sh
python3 local_runner/03_collect_outputs.py
```

Or run all at once:

```bash
bash local_runner/04_run_everything.sh
```

## Requirements

- Ubuntu/WSL shell
- `Rscript` available on PATH
- Python 3 for collector script
- Seurat object file present at one of:
  - `data/seu.qs`
  - `data/seu.rds`
  - `data/02_processed/placenta_infection_seurat.qs`
  - `data/02_processed/placenta_infection_seurat_scored.qs`
  - `seu_rna.full.rds`

## Output layout

Collected runs are written to:

- `analysis_runs/<timestamp>/figures`
- `analysis_runs/<timestamp>/tables`
- `analysis_runs/<timestamp>/logs`
- `analysis_runs/<timestamp>/objects`
- `analysis_runs/<timestamp>/manifest.tsv`

This does **not** delete your original `outputs/` or `results/`; it only organizes a clean copy.
