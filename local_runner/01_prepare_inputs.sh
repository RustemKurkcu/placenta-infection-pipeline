#!/usr/bin/env bash
set -euo pipefail

mkdir -p data

if [[ -f data/seu.qs || -f data/seu.rds ]]; then
  echo "[ok] data/seu.qs or data/seu.rds already present"
else
  if [[ -f data/02_processed/placenta_infection_seurat.qs ]]; then
    ln -sf ../data/02_processed/placenta_infection_seurat.qs data/seu.qs
    echo "[ok] linked data/seu.qs -> data/02_processed/placenta_infection_seurat.qs"
  elif [[ -f data/02_processed/placenta_infection_seurat_scored.qs ]]; then
    ln -sf ../data/02_processed/placenta_infection_seurat_scored.qs data/seu.qs
    echo "[ok] linked data/seu.qs -> data/02_processed/placenta_infection_seurat_scored.qs"
  elif [[ -f seu_rna.full.rds ]]; then
    ln -sf ../seu_rna.full.rds data/seu.rds
    echo "[ok] linked data/seu.rds -> seu_rna.full.rds"
  else
    echo "[error] no usable seurat object found" >&2
    exit 1
  fi
fi

for f in gene_sets/glyco_genes.txt gene_sets/adhesion_genes.txt gene_sets/innate_genes.txt; do
  if [[ ! -s "$f" ]]; then
    echo "[error] missing or empty required gene-set file: $f" >&2
    exit 1
  fi
done

echo "[ok] input preparation complete"
