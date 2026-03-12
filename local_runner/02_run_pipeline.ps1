$ErrorActionPreference = 'Stop'

$scripts = @(
    'scripts/01_load_make_seurat.R',
    'scripts/02_qc_overview.R',
    'scripts/03_core_umaps_and_composition.R',
    'scripts/04_gene_sets_scores_plots.R',
    'scripts/05_susceptibility_severity_models.R',
    'scripts/06_nk_cytotoxic_module.R',
    'scripts/07_integration_sensitivity.R',
    'scripts/08_organoid_vs_placenta_comparison.R',
    'scripts/09_reproducibility_report.R'
)

New-Item -ItemType Directory -Path 'outputs/logs' -Force | Out-Null
$masterLog = 'outputs/logs/local_runner_master.log'

$rscriptCmd = Get-Command 'Rscript.exe' -ErrorAction SilentlyContinue
if ($rscriptCmd) {
    $rscript = $rscriptCmd.Source
} else {
    $found = Get-ChildItem -Path 'C:\Program Files\R' -Filter 'Rscript.exe' -Recurse -ErrorAction SilentlyContinue | Select-Object -First 1
    if ($found) {
        $rscript = $found.FullName
    } else {
        throw '[error] Rscript.exe not found. Install R and ensure Rscript.exe is on PATH (or in C:\Program Files\R\...).'
    }
}

"[$(Get-Date -Format 's')] starting canonical pipeline" | Tee-Object -FilePath $masterLog -Append | Out-Host

foreach ($s in $scripts) {
    if (-not (Test-Path $s)) { throw "[error] missing script: $s" }

    $base = [System.IO.Path]::GetFileNameWithoutExtension($s)
    $stepLog = "outputs/logs/$base.log"

    "[$(Get-Date -Format 's')] RUN $s" | Tee-Object -FilePath $masterLog -Append | Out-Host

    $cmd = '"{0}" "{1}" > "{2}" 2>&1' -f $rscript, $s, $stepLog
    cmd.exe /c $cmd | Out-Null
    $exitCode = $LASTEXITCODE

    if ($exitCode -eq 0) {
        "[$(Get-Date -Format 's')] OK  $s" | Tee-Object -FilePath $masterLog -Append | Out-Host
    } else {
        "[$(Get-Date -Format 's')] FAIL $s (exit $exitCode; see $stepLog)" | Tee-Object -FilePath $masterLog -Append | Out-Host
        exit $exitCode
    }
}

"[$(Get-Date -Format 's')] pipeline complete" | Tee-Object -FilePath $masterLog -Append | Out-Host
