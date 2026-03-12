$ErrorActionPreference = 'Stop'

New-Item -ItemType Directory -Path 'data' -Force | Out-Null

function Ensure-LinkOrCopy {
    param(
        [Parameter(Mandatory=$true)][string]$Source,
        [Parameter(Mandatory=$true)][string]$Dest
    )

    if (Test-Path $Dest) { Remove-Item $Dest -Force }

    try {
        New-Item -ItemType SymbolicLink -Path $Dest -Target $Source -Force | Out-Null
        Write-Host "[ok] linked $Dest -> $Source"
    } catch {
        Copy-Item -Path $Source -Destination $Dest -Force
        Write-Host "[ok] copied $Source -> $Dest (symlink unavailable)"
    }
}

if ((Test-Path 'data/seu.qs') -or (Test-Path 'data/seu.rds')) {
    Write-Host '[ok] data/seu.qs or data/seu.rds already present'
} else {
    if (Test-Path 'data/02_processed/placenta_infection_seurat.qs') {
        Ensure-LinkOrCopy -Source 'data/02_processed/placenta_infection_seurat.qs' -Dest 'data/seu.qs'
    } elseif (Test-Path 'data/02_processed/placenta_infection_seurat_scored.qs') {
        Ensure-LinkOrCopy -Source 'data/02_processed/placenta_infection_seurat_scored.qs' -Dest 'data/seu.qs'
    } elseif (Test-Path 'seu_rna.full.rds') {
        Ensure-LinkOrCopy -Source 'seu_rna.full.rds' -Dest 'data/seu.rds'
    } else {
        throw '[error] no usable seurat object found'
    }
}

$required = @(
    'gene_sets/glyco_genes.txt',
    'gene_sets/adhesion_genes.txt',
    'gene_sets/innate_genes.txt'
)

foreach ($f in $required) {
    if (-not (Test-Path $f)) { throw "[error] missing required gene-set file: $f" }
    $fi = Get-Item $f
    if ($fi.Length -le 0) { throw "[error] empty required gene-set file: $f" }
}

Write-Host '[ok] input preparation complete'
