#!/usr/bin/env bash
set -euo pipefail

if ! command -v Rscript >/dev/null 2>&1; then
  echo "[error] Rscript not found on PATH" >&2
  exit 1
fi

mkdir -p outputs/logs
master_log="outputs/logs/local_runner_master.log"

scripts=(
  scripts/01_load_make_seurat.R
  scripts/02_qc_overview.R
  scripts/03_core_umaps_and_composition.R
  scripts/04_gene_sets_scores_plots.R
  scripts/05_susceptibility_severity_models.R
  scripts/06_nk_cytotoxic_module.R
  scripts/07_integration_sensitivity.R
  scripts/08_organoid_vs_placenta_comparison.R
  scripts/09_reproducibility_report.R
)

echo "[$(date '+%F %T')] starting canonical pipeline" | tee -a "$master_log"

for s in "${scripts[@]}"; do
  if [[ ! -f "$s" ]]; then
    echo "[$(date '+%F %T')] missing script: $s" | tee -a "$master_log"
    exit 1
  fi
  step_log="outputs/logs/$(basename "$s" .R).log"
  echo "[$(date '+%F %T')] RUN $s" | tee -a "$master_log"
  if Rscript "$s" >"$step_log" 2>&1; then
    echo "[$(date '+%F %T')] OK  $s" | tee -a "$master_log"
  else
    echo "[$(date '+%F %T')] FAIL $s (see $step_log)" | tee -a "$master_log"
    exit 1
  fi
done

echo "[$(date '+%F %T')] pipeline complete" | tee -a "$master_log"
