$ErrorActionPreference = 'Stop'

Write-Host '1) Preparing inputs...'
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\local_runner\01_prepare_inputs.ps1
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host '2) Running R pipeline...'
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\local_runner\02_run_pipeline.ps1
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host '3) Collecting outputs...'
python .\local_runner\03_collect_outputs.py
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host '[ok] full local run complete'
