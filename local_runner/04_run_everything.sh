#!/usr/bin/env bash
set -euo pipefail

bash local_runner/01_prepare_inputs.sh
bash local_runner/02_run_pipeline.sh
python3 local_runner/03_collect_outputs.py

echo "[ok] full local run complete"
