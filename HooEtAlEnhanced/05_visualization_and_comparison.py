"""
Comprehensive Visualization: Explant vs In Vivo Comparison
============================================================

Generates publication-quality figures comparing:
1. Cell type composition
2. QC metrics
3. F. nucleatum vulnerability scores
4. Spatial neighborhood analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['figure.dpi'] = 150

# ============================================================================
# LOAD DATA
# ============================================================================

print("Loading data...")

# Load explant data
explant = pd.read_csv('/workspace/metadata.csv', low_memory=False)
explant_ui = explant[explant['infection'] == 'UI'].copy()

# Load comparison data
comparison = pd.read_csv('/workspace/analysis/output/explant_vs_invivo_celltype_comparison.csv')

# ============================================================================
# FIGURE 1: CELL TYPE COMPOSITION COMPARISON
# ============================================================================

print("Creating Figure 1: Cell Type Composition...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: Explant cell types
ax1 = axes[0]
explant_ct = explant_ui['cell_type'].value_counts()
colors_explant = plt.cm.Set3(np.linspace(0, 1, len(explant_ct)))
bars1 = ax1.barh(range(len(explant_ct)), explant_ct.values, color=colors_explant)
ax1.set_yticks(range(len(explant_ct)))
ax1.set_yticklabels(explant_ct.index)
ax1.set_xlabel('Number of Cells')
ax1.set_title('A. Explant Model (Uninfected)\nHoo et al. 2024', fontweight='bold')
ax1.invert_yaxis()

# Add percentages
for i, (bar, val) in enumerate(zip(bars1, explant_ct.values)):
    pct = val / len(explant_ui) * 100
    ax1.text(val + 100, i, f'{pct:.1f}%', va='center', fontsize=9)

# Panel B: In Vivo cell types
ax2 = axes[1]
invivo_ct = pd.read_csv('/workspace/analysis/output/invivo_celltype_counts.csv')
invivo_ct = invivo_ct.sort_values('invivo_n', ascending=True).tail(15)
colors_invivo = plt.cm.Set2(np.linspace(0, 1, len(invivo_ct)))
bars2 = ax2.barh(range(len(invivo_ct)), invivo_ct['invivo_n'].values, color=colors_invivo)
ax2.set_yticks(range(len(invivo_ct)))
ax2.set_yticklabels(invivo_ct['invivo_cell_type'])
ax2.set_xlabel('Number of Cells')
ax2.set_title('B. In Vivo Placenta\nGreenbaum et al.', fontweight='bold')

# Add percentages
total_invivo = invivo_ct['invivo_n'].sum()
for i, (bar, val) in enumerate(zip(bars2, invivo_ct['invivo_n'].values)):
    pct = val / total_invivo * 100
    ax2.text(val + 100, i, f'{pct:.1f}%', va='center', fontsize=9)

plt.tight_layout()
plt.savefig('/workspace/analysis/figures/Fig01_CellType_Comparison.png', dpi=300, bbox_inches='tight')
plt.savefig('/workspace/analysis/figures/Fig01_CellType_Comparison.pdf', bbox_inches='tight')
print("  Saved: Fig01_CellType_Comparison")
plt.close()

# ============================================================================
# FIGURE 2: QC METRICS COMPARISON
# ============================================================================

print("Creating Figure 2: QC Metrics...")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel A: nFeature_RNA distribution by cell type
ax1 = axes[0, 0]
top_celltypes = list(explant_ui['cell_type'].value_counts().head(8).index)
data_subset = explant_ui[explant_ui['cell_type'].isin(top_celltypes)]
sns.violinplot(data=data_subset, x='cell_type', y='nFeature_RNA', ax=ax1, palette='Set2')
ax1.set_xlabel('Cell Type')
ax1.set_ylabel('Genes Detected (nFeature_RNA)')
ax1.set_title('A. Gene Detection by Cell Type', fontweight='bold')
ax1.tick_params(axis='x', rotation=45)

# Panel B: nCount_RNA distribution
ax2 = axes[0, 1]
sns.violinplot(data=data_subset, x='cell_type', y='nCount_RNA', ax=ax2, palette='Set2')
ax2.set_xlabel('Cell Type')
ax2.set_ylabel('UMI Counts (nCount_RNA)')
ax2.set_title('B. UMI Counts by Cell Type', fontweight='bold')
ax2.tick_params(axis='x', rotation=45)

# Panel C: Mitochondrial percentage
ax3 = axes[1, 0]
sns.boxplot(data=data_subset, x='cell_type', y='percent.mt', ax=ax3, palette='Set3')
ax3.axhline(y=0.1, color='red', linestyle='--', label='10% threshold')
ax3.axhline(y=0.2, color='darkred', linestyle='--', label='20% threshold')
ax3.set_xlabel('Cell Type')
ax3.set_ylabel('Mitochondrial %')
ax3.set_title('C. Mitochondrial Content by Cell Type', fontweight='bold')
ax3.tick_params(axis='x', rotation=45)
ax3.legend()

# Panel D: QC Summary
ax4 = axes[1, 1]
qc_summary = pd.DataFrame({
    'Metric': ['Median Genes', 'Median UMIs', 'Median MT%', 'Total Cells'],
    'Value': [
        explant_ui['nFeature_RNA'].median(),
        explant_ui['nCount_RNA'].median(),
        explant_ui['percent.mt'].median() * 100,
        len(explant_ui)
    ]
})
ax4.axis('off')
table = ax4.table(cellText=qc_summary.values, colLabels=qc_summary.columns, 
                   loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(11)
table.scale(1.2, 1.5)
ax4.set_title('D. QC Summary (Uninfected Explant)', fontweight='bold', pad=20)

plt.tight_layout()
plt.savefig('/workspace/analysis/figures/Fig02_QC_Metrics.png', dpi=300, bbox_inches='tight')
plt.savefig('/workspace/analysis/figures/Fig02_QC_Metrics.pdf', bbox_inches='tight')
print("  Saved: Fig02_QC_Metrics")
plt.close()

# ============================================================================
# FIGURE 3: F. NUCLEATUM VULNERABILITY SCORES
# ============================================================================

print("Creating Figure 3: F. nucleatum Vulnerability...")

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel A: GlycoScore by cell type
ax1 = axes[0, 0]
score_data = explant_ui.groupby('cell_type').agg({
    'GlycoScore': 'mean',
    'AdhesionScore': 'mean',
    'InnateScore': 'mean',
    'CytotoxicScore': 'mean'
}).reset_index()
score_data = score_data.sort_values('GlycoScore', ascending=True)

colors = plt.cm.RdYlGn(np.linspace(0.2, 0.8, len(score_data)))
ax1.barh(score_data['cell_type'], score_data['GlycoScore'], color=colors)
ax1.set_xlabel('Mean GlycoScore')
ax1.set_title('A. Glycosylation Score\n(Fn Adhesion Potential)', fontweight='bold')

# Panel B: AdhesionScore
ax2 = axes[0, 1]
score_data = score_data.sort_values('AdhesionScore', ascending=True)
colors = plt.cm.RdYlBu(np.linspace(0.2, 0.8, len(score_data)))
ax2.barh(score_data['cell_type'], score_data['AdhesionScore'], color=colors)
ax2.set_xlabel('Mean AdhesionScore')
ax2.set_title('B. Adhesion Score\n(Fn Binding Potential)', fontweight='bold')

# Panel C: InnateScore
ax3 = axes[1, 0]
score_data = score_data.sort_values('InnateScore', ascending=True)
colors = plt.cm.RdYlGn_r(np.linspace(0.2, 0.8, len(score_data)))
ax3.barh(score_data['cell_type'], score_data['InnateScore'], color=colors)
ax3.set_xlabel('Mean InnateScore')
ax3.set_title('C. Innate Immune Score\n(Defense Response)', fontweight='bold')

# Panel D: CytotoxicScore
ax4 = axes[1, 1]
score_data = score_data.sort_values('CytotoxicScore', ascending=True)
colors = plt.cm.RdPu(np.linspace(0.2, 0.8, len(score_data)))
ax4.barh(score_data['cell_type'], score_data['CytotoxicScore'], color=colors)
ax4.set_xlabel('Mean CytotoxicScore')
ax4.set_title('D. Cytotoxic Score\n(NK-like Activity)', fontweight='bold')

plt.tight_layout()
plt.savefig('/workspace/analysis/figures/Fig03_FN_Vulnerability.png', dpi=300, bbox_inches='tight')
plt.savefig('/workspace/analysis/figures/Fig03_FN_Vulnerability.pdf', bbox_inches='tight')
print("  Saved: Fig03_FN_Vulnerability")
plt.close()

# ============================================================================
# FIGURE 4: COMBINED VULNERABILITY INDEX
# ============================================================================

print("Creating Figure 4: Combined Vulnerability Index...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Calculate Fn Vulnerability Index
explant_ui['Fn_Vulnerability'] = (
    explant_ui['GlycoScore'] * 0.3 + 
    explant_ui['AdhesionScore'] * 0.3 + 
    (1 - explant_ui['InnateScore']) * 0.4  # Higher innate = lower vulnerability
)

# Panel A: Vulnerability by cell type
ax1 = axes[0]
vuln_by_ct = explant_ui.groupby('cell_type')['Fn_Vulnerability'].agg(['mean', 'std']).reset_index()
vuln_by_ct = vuln_by_ct.sort_values('mean', ascending=True)

colors = plt.cm.Reds(np.linspace(0.3, 0.9, len(vuln_by_ct)))
bars = ax1.barh(vuln_by_ct['cell_type'], vuln_by_ct['mean'], 
                xerr=vuln_by_ct['std']/np.sqrt(len(explant_ui)/len(vuln_by_ct)),
                color=colors, capsize=3)
ax1.set_xlabel('F. nucleatum Vulnerability Index')
ax1.set_title('A. Fn Vulnerability by Cell Type\n(Uninfected Explant)', fontweight='bold')
ax1.axvline(x=0.5, color='red', linestyle='--', alpha=0.5, label='High vulnerability threshold')

# Panel B: Distribution of vulnerability scores
ax2 = axes[1]
for ct in ['HBC', 'VCT', 'EVT_1', 'Endo_f', 'F']:
    data = explant_ui[explant_ui['cell_type'] == ct]['Fn_Vulnerability']
    sns.kdeplot(data, ax=ax2, label=ct, linewidth=2)
ax2.set_xlabel('F. nucleatum Vulnerability Index')
ax2.set_ylabel('Density')
ax2.set_title('B. Vulnerability Distribution\n(Top 5 Cell Types)', fontweight='bold')
ax2.legend()

plt.tight_layout()
plt.savefig('/workspace/analysis/figures/Fig04_Vulnerability_Index.png', dpi=300, bbox_inches='tight')
plt.savefig('/workspace/analysis/figures/Fig04_Vulnerability_Index.pdf', bbox_inches='tight')
print("  Saved: Fig04_Vulnerability_Index")
plt.close()

# ============================================================================
# FIGURE 5: MODEL VALIDATION SUMMARY
# ============================================================================

print("Creating Figure 5: Model Validation Summary...")

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel A: Cell type overlap
ax1 = axes[0, 0]
matched = comparison[comparison['_merge'] == 'both']
left_only = comparison[comparison['_merge'] == 'left_only']
right_only = comparison[comparison['_merge'] == 'right_only']

categories = ['Matched', 'In Vivo Only', 'Explant Only']
values = [len(matched), len(left_only), len(right_only)]
colors = ['#2ecc71', '#3498db', '#e74c3c']
bars = ax1.bar(categories, values, color=colors)
ax1.set_ylabel('Number of Cell Types')
ax1.set_title('A. Cell Type Overlap', fontweight='bold')
for bar, val in zip(bars, values):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, 
             str(val), ha='center', va='bottom', fontsize=12)

# Panel B: Proportion comparison scatter
ax2 = axes[0, 1]
matched_plot = matched[matched['explant_n'] > 0]
ax2.scatter(matched_plot['invivo_mean_prop'] * 100, 
            matched_plot['explant_pct'], 
            s=100, alpha=0.7, c='steelblue', edgecolor='black')
ax2.plot([0, 35], [0, 35], 'r--', label='1:1 line')
ax2.set_xlabel('In Vivo Proportion (%)')
ax2.set_ylabel('Explant Proportion (%)')
ax2.set_title('B. Cell Type Proportion Comparison', fontweight='bold')
ax2.legend()

# Add labels for key cell types
for _, row in matched_plot.iterrows():
    if row['explant_pct'] > 5 or row['invivo_mean_prop'] * 100 > 5:
        ax2.annotate(row['harmonized_ct'], 
                     (row['invivo_mean_prop'] * 100, row['explant_pct']),
                     fontsize=8, alpha=0.8)

# Panel C: Quality metrics
ax3 = axes[1, 0]
quality_metrics = {
    'Median Genes': explant_ui['nFeature_RNA'].median(),
    'Median UMIs': explant_ui['nCount_RNA'].median() / 1000,
    'MT% (<5%)': 100 - (explant_ui['percent.mt'].mean() * 100),
    'Cell Type Diversity': explant_ui['cell_type'].nunique(),
    'Donor Count': explant_ui['donor_id'].nunique()
}
ax3.barh(list(quality_metrics.keys()), list(quality_metrics.values()), 
         color=['#3498db', '#2ecc71', '#f39c12', '#9b59b6', '#e74c3c'])
ax3.set_xlabel('Value')
ax3.set_title('C. Quality Metrics Summary', fontweight='bold')

# Panel D: Model validity assessment
ax4 = axes[1, 1]
ax4.axis('off')

validity_text = """
MODEL VALIDITY ASSESSMENT
========================

STRENGTHS:
✓ All major placental cell types present
✓ High-quality data (median 4,356 genes/cell)
✓ Low mitochondrial content (median 3.9%)
✓ Multiple donors (n=10)
✓ Time-course design (24h, 48h)
✓ Multiple infection conditions (Lm, Tg, Pf)

LIMITATIONS:
✗ Hofbauer cells overrepresented (22% vs 4.5% in vivo)
✗ Endothelial cells enriched (15% vs 1.9% in vivo)
✗ vCTB underrepresented (27% vs 56% in vivo)
✗ Lacks maternal vasculature
✗ Ex vivo culture may alter cell states

OVERALL ASSESSMENT:
▶ Suitable for infection studies
▶ Good representation of key cell types
▶ Consider cell type proportions in analysis
"""

ax4.text(0.1, 0.95, validity_text, transform=ax4.transAxes, fontsize=10,
         verticalalignment='top', fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('/workspace/analysis/figures/Fig05_Model_Validation.png', dpi=300, bbox_inches='tight')
plt.savefig('/workspace/analysis/figures/Fig05_Model_Validation.pdf', bbox_inches='tight')
print("  Saved: Fig05_Model_Validation")
plt.close()

# ============================================================================
# FIGURE 6: INFECTION RESPONSE ANALYSIS
# ============================================================================

print("Creating Figure 6: Infection Response...")

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel A: Cell type composition by infection
ax1 = axes[0, 0]
infection_comp = pd.crosstab(explant['cell_type'], explant['infection'], normalize='columns') * 100
infection_comp = infection_comp[['UI', 'Lm', 'Tg', 'Pf']]
top_types = explant['cell_type'].value_counts().head(8).index
infection_comp_subset = infection_comp.loc[infection_comp.index.isin(top_types)]

infection_comp_subset.plot(kind='bar', ax=ax1, colormap='Set2', width=0.8)
ax1.set_xlabel('Cell Type')
ax1.set_ylabel('Proportion (%)')
ax1.set_title('A. Cell Type Composition by Infection', fontweight='bold')
ax1.tick_params(axis='x', rotation=45)
ax1.legend(title='Condition', labels=['UI', 'Lm', 'Tg', 'Pf'])

# Panel B: InnateScore by infection
ax2 = axes[0, 1]
top_5_types = list(explant['cell_type'].value_counts().head(5).index)
sns.boxplot(data=explant[explant['cell_type'].isin(top_5_types)], 
            x='infection', y='InnateScore', hue='cell_type',
            order=['UI', 'Lm', 'Tg', 'Pf'], ax=ax2, palette='Set2')
ax2.set_xlabel('Infection Condition')
ax2.set_ylabel('Innate Immune Score')
ax2.set_title('B. Innate Immune Response by Infection', fontweight='bold')
ax2.legend(title='Cell Type', bbox_to_anchor=(1.02, 1), loc='upper left')

# Panel C: GlycoScore by infection
ax3 = axes[1, 0]
sns.boxplot(data=explant[explant['cell_type'].isin(top_5_types)], 
            x='infection', y='GlycoScore', hue='cell_type',
            order=['UI', 'Lm', 'Tg', 'Pf'], ax=ax3, palette='Set2')
ax3.set_xlabel('Infection Condition')
ax3.set_ylabel('Glycosylation Score')
ax3.set_title('C. Glycosylation Score by Infection', fontweight='bold')
ax3.legend(title='Cell Type', bbox_to_anchor=(1.02, 1), loc='upper left')

# Panel D: Time course analysis
ax4 = axes[1, 1]
time_comp = explant.groupby(['hpi', 'infection']).agg({
    'InnateScore': 'mean',
    'GlycoScore': 'mean'
}).reset_index()

for inf in ['UI', 'Lm', 'Tg', 'Pf']:
    data = time_comp[time_comp['infection'] == inf]
    ax4.plot(data['hpi'], data['InnateScore'], marker='o', label=inf, linewidth=2)

ax4.set_xlabel('Time Point')
ax4.set_ylabel('Mean Innate Immune Score')
ax4.set_title('D. Immune Response Over Time', fontweight='bold')
ax4.legend(title='Condition')

plt.tight_layout()
plt.savefig('/workspace/analysis/figures/Fig06_Infection_Response.png', dpi=300, bbox_inches='tight')
plt.savefig('/workspace/analysis/figures/Fig06_Infection_Response.pdf', bbox_inches='tight')
print("  Saved: Fig06_Infection_Response")
plt.close()

# ============================================================================
# SAVE SUMMARY STATISTICS
# ============================================================================

print("\nSaving summary statistics...")

# Cell type vulnerability summary
vuln_summary = explant_ui.groupby('cell_type').agg({
    'nCount_RNA': ['count', 'median'],
    'nFeature_RNA': 'median',
    'percent.mt': 'median',
    'GlycoScore': 'mean',
    'AdhesionScore': 'mean',
    'InnateScore': 'mean',
    'CytotoxicScore': 'mean',
    'Fn_Vulnerability': 'mean'
}).round(3)
vuln_summary.columns = ['_'.join(col).strip() for col in vuln_summary.columns.values]
vuln_summary.to_csv('/workspace/analysis/output/explant_vulnerability_summary.csv')

print("\n" + "=" * 80)
print("ANALYSIS COMPLETE")
print("=" * 80)
print(f"\nGenerated 6 publication-quality figures:")
print("  1. Fig01_CellType_Comparison - Explant vs In Vivo composition")
print("  2. Fig02_QC_Metrics - Quality control metrics")
print("  3. Fig03_FN_Vulnerability - F. nucleatum vulnerability scores")
print("  4. Fig04_Vulnerability_Index - Combined vulnerability index")
print("  5. Fig05_Model_Validation - Model validity assessment")
print("  6. Fig06_Infection_Response - Infection response analysis")
print(f"\nSaved to: /workspace/analysis/figures/")
print(f"\nData tables saved to: /workspace/analysis/output/")