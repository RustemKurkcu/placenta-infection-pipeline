#!/usr/bin/env python3
"""
01_prepare_invivo_data.py
=========================
Prepare BOTH in vivo placenta datasets from Greenbaum et al. (hPlacenta-architecture):
1. Slide-tags spatial (1,923 cells, 10 cell types) - has expression + spatial coords
2. Multiome reference (36,456 cells, 17 cell types) - full reference (metadata only here)

For comparison with Hoo et al. 2024 organoid/explant data.
"""

import pandas as pd
import numpy as np
import gzip
import os

DATA_DIR = "hPlacenta-architecture/data/raw/Broad_SCP2601human-placenta-architecture"
OUTPUT_DIR = "analysis/output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 70)
print("PART A: Slide-tags Spatial Data (1,923 cells with expression)")
print("=" * 70)

# Load cluster annotations (Slide-tags)
cluster_path = os.path.join(DATA_DIR, "humanplacenta_cluster.csv")
cluster_df = pd.read_csv(cluster_path, skiprows=[1])  # Skip TYPE row
cluster_df.columns = [c.strip('"') for c in cluster_df.columns]
cluster_df["NAME"] = cluster_df["NAME"].str.strip('"')
cluster_df["cell_type"] = cluster_df["cell_type"].str.strip('"')
cluster_df = cluster_df.set_index("NAME")
print(f"Cluster annotations: {len(cluster_df)} cells")

# Load spatial coordinates
spatial_path = os.path.join(DATA_DIR, "humanplacenta_spatial.csv")
spatial_df = pd.read_csv(spatial_path, skiprows=[1])
spatial_df.columns = [c.strip('"') for c in spatial_df.columns]
spatial_df["NAME"] = spatial_df["NAME"].str.strip('"')
spatial_df = spatial_df.set_index("NAME")
# Convert X, Y to numeric
for col in ["X", "Y"]:
    if spatial_df[col].dtype == object:
        spatial_df[col] = pd.to_numeric(spatial_df[col].str.strip('"'), errors='coerce')
    else:
        spatial_df[col] = pd.to_numeric(spatial_df[col], errors='coerce')
print(f"Spatial coordinates: {len(spatial_df)} cells")

# Cell type distribution
print("\n--- Slide-tags Cell Type Distribution ---")
ct_counts = cluster_df["cell_type"].value_counts()
for ct, count in ct_counts.items():
    print(f"  {ct:25s}: {count:5d} ({100*count/len(cluster_df):.1f}%)")

# Load expression data
print("\nLoading expression matrix...")
expr_path = os.path.join(DATA_DIR, "humanplacenta_expression.csv.gz")
expr = pd.read_csv(expr_path, index_col=0)
print(f"Expression: {expr.shape[0]} genes x {expr.shape[1]} cells")

# Verify overlap
common = set(expr.columns) & set(cluster_df.index)
print(f"Cells in both: {len(common)}")

# Compute per-cell-type mean expression
print("\n--- Computing Per-Cell-Type Mean Expression ---")
celltype_groups = cluster_df.loc[cluster_df.index.isin(common)].groupby("cell_type").groups

mean_expr = {}
detect_rate = {}
cell_counts = {}
for ct, cells in celltype_groups.items():
    cells_list = [c for c in cells if c in expr.columns]
    if len(cells_list) >= 5:  # minimum 5 cells
        sub = expr[cells_list]
        mean_expr[ct] = sub.mean(axis=1)
        detect_rate[ct] = (sub > 0).mean(axis=1)
        cell_counts[ct] = len(cells_list)
        print(f"  {ct:25s}: {len(cells_list):5d} cells")

mean_expr_df = pd.DataFrame(mean_expr)
detect_rate_df = pd.DataFrame(detect_rate)
print(f"\nMean expression: {mean_expr_df.shape[0]} genes x {mean_expr_df.shape[1]} cell types")

# Save
mean_expr_df.to_csv(os.path.join(OUTPUT_DIR, "slidetags_mean_expression_by_celltype.csv"))
detect_rate_df.to_csv(os.path.join(OUTPUT_DIR, "slidetags_detection_rate_by_celltype.csv"))

print("\n" + "=" * 70)
print("PART B: Multiome Reference Metadata (36,456 cells)")
print("=" * 70)

meta_path = os.path.join(DATA_DIR, "metadata.csv")
meta = pd.read_csv(meta_path, skiprows=[1])
meta = meta.set_index("NAME")
print(f"Multiome metadata: {len(meta)} cells")

ct_counts_multi = meta["cluster"].value_counts()
print("\n--- Multiome Cell Type Distribution ---")
for ct, count in ct_counts_multi.items():
    print(f"  {ct:30s}: {count:6d} ({100*count/len(meta):.1f}%)")

# Donor distribution
print("\n--- Donor Distribution ---")
for d, c in meta["biosample_id"].value_counts().items():
    print(f"  {d}: {c} cells")

# Save multiome summary
meta_summary = pd.DataFrame({
    "cell_type": ct_counts_multi.index,
    "count": ct_counts_multi.values,
    "pct": (100 * ct_counts_multi.values / len(meta)).round(1)
})
meta_summary.to_csv(os.path.join(OUTPUT_DIR, "multiome_celltype_distribution.csv"), index=False)

print("\n" + "=" * 70)
print("PART C: Key Gene Expression Analysis")
print("=" * 70)

# Define key gene sets
gene_sets = {
    "MMP_ECM_Remodeling": ["MMP2","MMP9","MMP14","ADAMTS4","ADAMTS5","CTSK","CTSL","CTSB",
                           "TIMP1","TIMP2","SERPINE1","PLAUR"],
    "ECM_Structure": ["COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","COL6A2","FN1","DCN","LUM"],
    "Immune_Tolerance": ["HLA-G","CD274","PDCD1LG2","IDO1","TGFB1","IL10","LGALS9","HAVCR2",
                         "VSIR","LILRB1","LILRB2"],
    "Cytotoxic_NK": ["NKG7","GNLY","PRF1","GZMB","GZMA","GZMH","CTSW","XCL1","XCL2",
                     "KLRD1","KLRK1","NCR1","FCGR3A","TYROBP","TRAC"],
    "Ethanolamine_Metabolism": ["PLD1","ETNK1","ETNK2","PCYT2","SELENOI","CEPT1","PLA2G6"],
    "EA_Availability": ["GDPD1","GDPD5","PLD1","FAAH","ETNK1","PCYT2"],
    "Endothelial_Activation": ["VCAM1","ICAM1","SELE","ENG","EDN1","NOS3","KDR","FLT1","PGF","VEGFA"],
    "Hypoxia_Response": ["HIF1A","VEGFA","SLC2A1","LDHA","CA9","BNIP3","EGLN1"],
    "Entry_Receptors": ["CDH1","EPCAM","GALNT1"],
    "Inflammatory": ["IL1B","TNF","NFKBIA","CXCL8","S100A8","S100A9","CCL2","CCL3","CCL4"],
    "Interferon_Stimulated": ["ISG15","IFIT1","IFIT2","IFIT3","MX1","OAS1","OAS2","OAS3","IFI6","IFI27"],
    "Trophoblast_Core": ["KRT7","KRT8","KRT18","TACSTD2","PAGE4","GCM1","CGA"],
    "EVT_Specific": ["HLA-G","ITGA1","MMP2","MMP9","ITGA5","NOTCH1","NOTCH2"],
    "HBC_Macrophage": ["CD68","CD163","CSF1R","C1QA","C1QB","C1QC","APOE","TYROBP","LST1"],
    "Fibroblast_Core": ["COL1A1","COL1A2","DCN","LUM","COL3A1","VIM","ACTA2"],
    "Endothelial_Core": ["PECAM1","VWF","KDR","ESAM","RAMP2","CDH5"],
    "Glycosylation_Fn": ["GALNT1","GALNT2","GALNT3","GALNT6","GALNT7","GALNT10",
                         "B3GNT5","B4GALT1","FUT2","FUT4","FUT8","ST3GAL1","ST6GAL1"],
}

all_genes = set(mean_expr_df.index)
print(f"\nTotal genes in Slide-tags expression: {len(all_genes)}")

print("\nGene Set Availability:")
for gs_name, genes in gene_sets.items():
    found = [g for g in genes if g in all_genes]
    missing = [g for g in genes if g not in all_genes]
    pct = 100 * len(found) / len(genes) if genes else 0
    status = "✓" if pct == 100 else ("◐" if pct >= 50 else "✗")
    print(f"  {status} {gs_name:30s}: {len(found):2d}/{len(genes):2d} ({pct:.0f}%)", end="")
    if missing and len(missing) <= 5:
        print(f"  missing: {', '.join(missing)}")
    elif missing:
        print(f"  missing: {', '.join(missing[:3])}...")
    else:
        print()

# Extract key gene expression profiles
all_key_genes = sorted(set(g for genes in gene_sets.values() for g in genes))
found_key = [g for g in all_key_genes if g in all_genes]
print(f"\nKey genes found: {len(found_key)}/{len(all_key_genes)}")

key_expr = mean_expr_df.loc[found_key]
key_detect = detect_rate_df.loc[found_key]
key_expr.to_csv(os.path.join(OUTPUT_DIR, "slidetags_key_genes_mean_expr.csv"))
key_detect.to_csv(os.path.join(OUTPUT_DIR, "slidetags_key_genes_detect_rate.csv"))

# EA Availability Score
print("\n--- EA Availability Score per Cell Type ---")
eas_num = [g for g in ["GDPD1","GDPD5","PLD1","FAAH"] if g in all_genes]
eas_den = [g for g in ["ETNK1","PCYT2"] if g in all_genes]
print(f"  Numerator genes found: {eas_num}")
print(f"  Denominator genes found: {eas_den}")

if eas_num and eas_den:
    for ct in mean_expr_df.columns:
        num = mean_expr_df.loc[eas_num, ct].sum()
        den = mean_expr_df.loc[eas_den, ct].sum()
        eas = num / (den + 0.01)
        print(f"    {ct:25s}: EAS = {eas:.4f} (num={num:.4f}, den={den:.4f})")

# Top markers per cell type
print("\n--- Top 15 Marker Genes per Cell Type (by mean expression) ---")
for ct in mean_expr_df.columns:
    top = mean_expr_df[ct].nlargest(15)
    genes_str = ", ".join([f"{g}({v:.1f})" for g, v in top.items()])
    print(f"\n  {ct}:")
    print(f"    {genes_str}")

# Gene variance across cell types
gene_var = mean_expr_df.var(axis=1).sort_values(ascending=False)
print("\n--- Top 30 Most Variable Genes Across Cell Types ---")
for g, v in gene_var.head(30).items():
    # Which cell type has max expression?
    max_ct = mean_expr_df.loc[g].idxmax()
    max_val = mean_expr_df.loc[g].max()
    print(f"  {g:20s}: var={v:.3f}, max in {max_ct} ({max_val:.2f})")

gene_var.to_csv(os.path.join(OUTPUT_DIR, "slidetags_gene_variance.csv"))

print("\n" + "=" * 70)
print("PART D: Cell Type Mapping for Organoid Comparison")
print("=" * 70)

# Map between Greenbaum (in vivo) and Hoo et al. (organoid/explant) cell types
celltype_mapping = {
    # Greenbaum Slide-tags -> Hoo et al. explant
    "vCTB": ["VCT", "VCT_fusing", "VCT_p", "VCT_CCC"],
    "STB": ["SCT"],
    "STB-progenitor": ["VCT_fusing"],
    "EVT": ["iEVT", "EVT_1", "EVT_2"],
    "EVT-progenitor": ["iEVT"],
    "FIB1": ["F", "F_p", "F_sm"],
    "FIB2": ["F", "F_p"],
    "Hofbauer cells": ["HBC", "HBC_p"],
    "Endothelial": ["Endo_f"],
    "Erythroblasts": [],  # No direct match
}

# Also map multiome cell types
multiome_mapping = {
    "vCTB1": ["VCT"],
    "vCTB2": ["VCT"],
    "vCTB3": ["VCT"],
    "vCTBp": ["VCT_p"],
    "STB": ["SCT"],
    "EVT1": ["EVT_1"],
    "EVT2": ["EVT_2"],
    "EVT3": ["iEVT"],
    "Fibroblast1": ["F"],
    "Fibroblast2": ["F"],
    "Endothelial_cells": ["Endo_f"],
    "Hofbauer cells": ["HBC"],
    "maternal macrophages ": ["PAMM1"],
    "Maternal fibroblast": [],
}

print("\nCell Type Mapping (In Vivo → Organoid/Explant):")
print("\n  Slide-tags (spatial):")
for invivo, explant in celltype_mapping.items():
    if explant:
        print(f"    {invivo:25s} → {', '.join(explant)}")
    else:
        print(f"    {invivo:25s} → (no direct match)")

print("\n  Multiome (reference):")
for invivo, explant in multiome_mapping.items():
    if explant:
        print(f"    {invivo:30s} → {', '.join(explant)}")

# Save mapping
mapping_records = []
for invivo, explant_list in celltype_mapping.items():
    for e in explant_list:
        mapping_records.append({"invivo_slidetags": invivo, "explant_hoo": e})
pd.DataFrame(mapping_records).to_csv(os.path.join(OUTPUT_DIR, "celltype_mapping.csv"), index=False)

print("\n✅ All in vivo data prepared successfully!")
print(f"Outputs saved to: {OUTPUT_DIR}/")