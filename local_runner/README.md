# Local Runner Bundle (manifest-aware)

This folder gives you a clean, reproducible way to run the pipeline on your own workstation and collect outputs in one logical place.

## What this bundle does

1. **Prepare inputs** from the files you already have (large `.qs`/`.rds` objects).
2. **Run canonical scripts in order** (`scripts/01..09`) with per-step logs.
3. **Collect figures/tables/logs** from both `outputs/` and `results/` into one timestamped run folder.

## Recommended run order

### Linux / WSL / Git-Bash

```bash
bash local_runner/01_prepare_inputs.sh
bash local_runner/02_run_pipeline.sh
python3 local_runner/03_collect_outputs.py
```

Or all at once:

```bash
bash local_runner/04_run_everything.sh
```

### Windows PowerShell (recommended on your setup)

```powershell
powershell.exe -ExecutionPolicy Bypass -File .\local_runner\01_prepare_inputs.ps1
powershell.exe -ExecutionPolicy Bypass -File .\local_runner\02_run_pipeline.ps1
python .\local_runner\03_collect_outputs.py
```

Or all at once:

```powershell
powershell.exe -ExecutionPolicy Bypass -File .\local_runner\04_run_everything.ps1
```

## Requirements

- `Rscript` available on PATH (or installed under `C:\Program Files\R\...`)
- Python 3 for collector script
- Seurat object file present at one of:
  - `data/seu.qs`
  - `data/seu.rds`
  - `data/02_processed/placenta_infection_seurat.qs`
  - `data/02_processed/placenta_infection_seurat_scored.qs`
  - `seu_rna.full.rds`

## Notes for the issue you hit

- The previous PowerShell one-shot script could continue to collection even when R failed; this is now fixed (`04_run_everything.ps1` exits on pipeline failure).
- `02_run_pipeline.ps1` now captures stdout+stderr via `cmd.exe` redirection, avoiding the `Start-Process` same-file redirection error.
- `Rscript.exe : Warning messages:` lines are not automatically fatal. Check the script-specific log in `outputs/logs/*.log` and the master status in `outputs/logs/local_runner_master.log`.

## Output layout

Collected runs are written to:

- `analysis_runs/<timestamp>/figures`
- `analysis_runs/<timestamp>/tables`
- `analysis_runs/<timestamp>/logs`
- `analysis_runs/<timestamp>/objects`
- `analysis_runs/<timestamp>/manifest.tsv`

This does **not** delete your original `outputs/` or `results/`; it only organizes a clean copy.
