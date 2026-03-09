# Git LFS history cleanup for oversized `.qs` files

If `git-filter-repo` fails with:

> `Error: Option --paths unrecognized; did you mean --path or --paths-from-file?`

use `--path` (singular) for each file.

## 1) Rewrite history to remove >2 GiB files

Run from repo root:

```powershell
git checkout main

# Use git-filter-repo if on PATH; otherwise python -m git_filter_repo
$gfr = (Get-Command git-filter-repo -ErrorAction SilentlyContinue).Source
if ($gfr) {
  & $gfr --invert-paths `
    --path "data/02_processed/placenta_infection_seurat_scored.qs" `
    --path "data/02_processed/placenta_infection_seurat.qs" `
    --path "data/seu_clean.qs" `
    --path "outputs/objects/seu_with_nk_flags.qs" `
    --path "outputs/objects/seu_with_scores.qs" `
    --path "outputs/objects/seu_core_figs.qs" `
    --path "outputs/objects/seu_qc_checked.qs"
} else {
  python -m git_filter_repo --invert-paths `
    --path "data/02_processed/placenta_infection_seurat_scored.qs" `
    --path "data/02_processed/placenta_infection_seurat.qs" `
    --path "data/seu_clean.qs" `
    --path "outputs/objects/seu_with_nk_flags.qs" `
    --path "outputs/objects/seu_with_scores.qs" `
    --path "outputs/objects/seu_core_figs.qs" `
    --path "outputs/objects/seu_qc_checked.qs"
}
```

## 2) Cleanup local object references

```powershell
git reflog expire --expire=now --all
git gc --prune=now --aggressive
git lfs prune
```

## 3) Verify + push rewritten history

```powershell
# Informational: files > 2 GiB in working tree
Get-ChildItem -Recurse -File | Where-Object { $_.Length -gt 2147483648 } |
  Format-Table FullName, @{Name='GB';Expression={[math]::Round($_.Length/1GB,3)}} -AutoSize

# LFS pointers still tracked
git lfs ls-files

# Publish rewritten history
git push origin --force --all
git push origin --force --tags
```

## 4) If push still fails

Capture and inspect objects still referenced by LFS:

```powershell
git lfs ls-files -l
```

Then remove remaining offending paths from history and rerun steps 2–3.
