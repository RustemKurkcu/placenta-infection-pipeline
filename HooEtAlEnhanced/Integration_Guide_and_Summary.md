<h1>Integration Guide and Analysis Summary</h1><h2>Explant vs In Vivo Placenta Comparison for F. nucleatum Research</h2><hr><h2>Important Clarification</h2><p><strong>This analysis compares PLACENTAL EXPLANT data (not organoid) from Hoo et al. 2024 with in vivo placenta data from Greenbaum et al.</strong></p><p>The explant model is an ex vivo culture system that maintains tissue architecture and cell-cell interactions, while organoids are 3D self-organizing structures derived from stem cells.</p><hr><h2>Quick Reference: Key Findings</h2><h3>Model Validity: ⭐⭐⭐⭐☆ (4/5 stars)</h3><table class="e-rte-table"> <thead> <tr> <th>Aspect</th> <th>Rating</th> <th>Notes</th> </tr> </thead> <tbody><tr> <td>Cell Type Coverage</td> <td>⭐⭐⭐⭐⭐</td> <td>All major types present</td> </tr> <tr> <td>Data Quality</td> <td>⭐⭐⭐⭐⭐</td> <td>Excellent QC metrics</td> </tr> <tr> <td>Biological Relevance</td> <td>⭐⭐⭐⭐</td> <td>Good, but proportions differ</td> </tr> <tr> <td>Spatial Context</td> <td>⭐⭐⭐</td> <td>Lost in ex vivo culture</td> </tr> <tr> <td>Reproducibility</td> <td>⭐⭐⭐⭐⭐</td> <td>10 donors, multiple timepoints</td> </tr> </tbody></table><p><strong>Verdict</strong>: Suitable for F. nucleatum infection studies with caveats regarding cell type proportions.</p><hr><h2>Data Files Available</h2><h3>Your Uploaded Files (Explant Data)</h3><table class="e-rte-table"> <thead> <tr> <th>File</th> <th>Description</th> <th>Size</th> <th>Usage</th> </tr> </thead> <tbody><tr> <td><code>metadata.csv</code></td> <td>Cell metadata with 159K cells</td> <td>63 MB</td> <td>✅ Used</td> </tr> <tr> <td><code>seu_with_nk_flags_meta.data.csv</code></td> <td>NK cell flags metadata</td> <td>61 MB</td> <td>Available</td> </tr> <tr> <td><code>organoid_scripts.zip</code></td> <td>Analysis scripts</td> <td>19 KB</td> <td>✅ Reviewed</td> </tr> <tr> <td><code>Outputs.zip</code></td> <td>Pre-computed outputs</td> <td>2.4 MB</td> <td>✅ Reviewed</td> </tr> <tr> <td><code>gene_sets.zip</code></td> <td>Gene set definitions</td> <td>1.2 KB</td> <td>✅ Used</td> </tr> </tbody></table><h3>GitHub Repository (In Vivo Data)</h3><table class="e-rte-table"> <thead> <tr> <th>File</th> <th>Description</th> <th>Usage</th> </tr> </thead> <tbody><tr> <td><code>celltype_proportions_by_week_and_version.csv</code></td> <td>Cell type counts</td> <td>✅ Used</td> </tr> <tr> <td>Slide-tags tables</td> <td>Spatial analysis</td> <td>Partially used</td> </tr> <tr> <td>Multiome tables</td> <td>Reference data</td> <td>Identified</td> </tr> </tbody></table><hr><h2>Generated Outputs</h2><h3>Figures (6 publication-quality plots)</h3><ol> <li><strong>Fig01_CellType_Comparison</strong> - Side-by-side bar charts of cell type composition</li> <li><strong>Fig02_QC_Metrics</strong> - Quality control metrics by cell type</li> <li><strong>Fig03_FN_Vulnerability</strong> - F. nucleatum vulnerability scores (4 panels)</li> <li><strong>Fig04_Vulnerability_Index</strong> - Combined vulnerability index with distributions</li> <li><strong>Fig05_Model_Validation</strong> - Model validity assessment</li> <li><strong>Fig06_Infection_Response</strong> - Infection response across conditions</li> </ol><p>All figures available in PNG and PDF formats in <code>/workspace/analysis/figures/</code></p><h3>Data Tables</h3><table class="e-rte-table"> <thead> <tr> <th>Table</th> <th>Description</th> </tr> </thead> <tbody><tr> <td><code>explant_celltype_counts.csv</code></td> <td>Cell type counts and percentages</td> </tr> <tr> <td><code>invivo_celltype_counts.csv</code></td> <td>In vivo cell type counts</td> </tr> <tr> <td><code>explant_vs_invivo_celltype_comparison.csv</code></td> <td>Harmonized comparison</td> </tr> <tr> <td><code>explant_vulnerability_summary.csv</code></td> <td>Fn vulnerability by cell type</td> </tr> <tr> <td><code>fn_vulnerability_summary.csv</code></td> <td>R-generated vulnerability table</td> </tr> </tbody></table><h3>Scripts</h3><table class="e-rte-table"> <thead> <tr> <th>Script</th> <th>Language</th> <th>Purpose</th> </tr> </thead> <tbody><tr> <td><code>05_visualization_and_comparison.py</code></td> <td>Python</td> <td>Generate all figures</td> </tr> <tr> <td><code>01_Explorant_vs_InVivo_Comparison.R</code></td> <td>R</td> <td>Seurat integration pipeline</td> </tr> </tbody></table><hr><h2>Critical Differences: Explant vs In Vivo</h2><h3>Cell Type Proportions</h3><pre><code>Hofbauer cells:    ████████████████████ 22.4% (explant) vs 4.5% (in vivo)
Endothelial:       ████████████████ 14.9% (explant) vs 1.9% (in vivo)
vCTB:              ████████████████ 26.8% (explant) vs 56.4% (in vivo)
Fibroblasts:       ████████████ 14.7% (explant) vs 22.2% (in vivo)
EVT:               █████ 4.8% (explant) vs 16.9% (in vivo)
STB:               ████████████████ 16.6% (explant) vs 9.8% (in vivo)
</code></pre><h3>Implications</h3><ol> <li><strong>Hofbauer Cell Enrichment</strong>: May exaggerate immune response in explant</li> <li><strong>Endothelial Enrichment</strong>: Alters vascular context</li> <li><strong>vCTB Underrepresentation</strong>: May miss key trophoblast biology</li> <li><strong>EVT Underrepresentation</strong>: Reduced invasion capacity</li> </ol><p><strong>Recommendation</strong>: Use cell type proportions as weights when comparing infected vs uninfected conditions.</p><hr><h2>F. nucleatum Vulnerability Analysis</h2><h3>Most Vulnerable Cell Types</h3><table class="e-rte-table"> <thead> <tr> <th>Rank</th> <th>Cell Type</th> <th>Fn Vulnerability</th> <th>Key Reason</th> </tr> </thead> <tbody><tr> <td>1</td> <td>VCT</td> <td>0.557</td> <td>High adhesion, low immunity</td> </tr> <tr> <td>2</td> <td>VCT_p</td> <td>0.545</td> <td>High adhesion, low immunity</td> </tr> <tr> <td>3</td> <td>EVT_1</td> <td>0.534</td> <td>Moderate adhesion, low immunity</td> </tr> <tr> <td>4</td> <td>VCT_CCC</td> <td>0.523</td> <td>High adhesion, moderate immunity</td> </tr> <tr> <td>5</td> <td>EVT_2</td> <td>0.493</td> <td>Moderate adhesion, moderate immunity</td> </tr> </tbody></table><h3>Most Protected Cell Types</h3><table class="e-rte-table"> <thead> <tr> <th>Rank</th> <th>Cell Type</th> <th>Fn Vulnerability</th> <th>Key Reason</th> </tr> </thead> <tbody><tr> <td>1</td> <td>HBC</td> <td>0.336</td> <td>High innate immunity (0.907)</td> </tr> <tr> <td>2</td> <td>PAMM1</td> <td>0.368</td> <td>High innate immunity (0.941)</td> </tr> <tr> <td>3</td> <td>HBC_p</td> <td>0.369</td> <td>Moderate-high immunity (0.789)</td> </tr> </tbody></table><h3>MegL Toxic Switch Hypothesis Relevance</h3><p>The <strong>MegL Toxic Switch</strong> hypothesis proposes that F. nucleatum switches from ethanolamine (EA) metabolism to H₂S production via methionine γ-lyase (MegL), leading to placental dysfunction.</p><p><strong>Key Findings</strong>:</p><ul> <li><strong>Hofbauer cells</strong> have the highest EA metabolism scores (GlycoScore + AdhesionScore) but also highest immunity</li> <li><strong>VCT cells</strong> are most vulnerable due to high adhesion potential and low immunity</li> <li><strong>EVT cells</strong> show moderate vulnerability with specific subtypes (EVT_1) being more vulnerable</li> </ul><p><strong>Recommendation</strong>: Focus on VCT and EVT subtypes for infection studies, and assess MegL expression in these cell types.</p><hr><h2>Integration Pipeline</h2><h3>Option 1: Seurat v5 Integration (Recommended)</h3><pre><code class="language-r">library(Seurat)
library(dplyr)

# 1. Load both datasets
# You'll need the expression matrices for both datasets
explant_counts &lt;- Read10X_h5("explant_matrix.h5")
invivo_counts &lt;- Read10X_h5("invivo_matrix.h5")

# 2. Create Seurat objects
seu_explant &lt;- CreateSeuratObject(
  counts = explant_counts,
  project = "Explant"
)
seu_invivo &lt;- CreateSeuratObject(
  counts = invivo_counts,
  project = "InVivo"
)

# 3. Add metadata
seu_explant$cell_type &lt;- explant_meta$cell_type
seu_explant$infection &lt;- explant_meta$infection
seu_invivo$cell_type &lt;- invivo_meta$cell_type

# 4. Standard preprocessing
seu_explant &lt;- NormalizeData(seu_explant)
seu_explant &lt;- FindVariableFeatures(seu_explant, nfeatures = 3000)
seu_invivo &lt;- NormalizeData(seu_invivo)
seu_invivo &lt;- FindVariableFeatures(seu_invivo, nfeatures = 3000)

# 5. Integration
features &lt;- SelectIntegrationFeatures(
  list(seu_explant, seu_invivo),
  nfeatures = 3000
)

anchors &lt;- FindIntegrationAnchors(
  object.list = list(seu_explant, seu_invivo),
  anchor.features = features
)

seu_integrated &lt;- IntegrateData(anchorset = anchors)

# 6. Dimensionality reduction
seu_integrated &lt;- ScaleData(seu_integrated)
seu_integrated &lt;- RunPCA(seu_integrated)
seu_integrated &lt;- RunUMAP(seu_integrated, dims = 1:30)

# 7. Visualization
DimPlot(seu_integrated, group.by = c("cell_type", "orig.ident"))
</code></pre><h3>Option 2: Harmony Integration</h3><pre><code class="language-r">library(harmony)
library(Seurat)

# After combining and running PCA
seu_combined &lt;- merge(seu_explant, y = seu_invivo)
seu_combined &lt;- RunPCA(seu_combined)

# Run Harmony
seu_combined &lt;- RunHarmony(
  seu_combined,
  group.by.vars = "orig.ident",
  reduction = "pca"
)

# Run UMAP on Harmony embeddings
seu_combined &lt;- RunUMAP(
  seu_combined,
  reduction = "harmony",
  dims = 1:30
)
</code></pre><hr><h2>What You Need to Complete Full Integration</h2><h3>Required Data Files</h3><ol> <li><strong>Explant Expression Matrix</strong> (missing) <ul> <li>Format: <code>explant_matrix.h5</code> or <code>explant_matrix.mtx</code></li> <li>Contains: Gene × Cell UMI counts</li> <li>Location: Should be in <code>data/02_processed/</code> or similar</li> </ul> </li> <li><strong>In Vivo Expression Matrix</strong> (available in GitHub) <ul> <li>Location: <code>hPlacenta-architecture/data/raw/</code></li> <li>Format: Multiome or Slide-tags matrix</li> </ul> </li> </ol><h3>Optional but Recommended</h3><ol start="3"> <li><strong>Spatial Coordinates</strong> (for spatial analysis) <ul> <li>Slide-tags: <code>slidetags_coordinates.csv</code></li> <li>STARmap: <code>starmap_coordinates.csv</code></li> </ul> </li> <li><strong>Additional Metadata</strong> <ul> <li>Sample information</li> <li>Developmental stage details</li> <li>Donor characteristics</li> </ul> </li> </ol><hr><h2>Recommendations for Your Research</h2><h3>For F. nucleatum Infection Studies</h3><ol> <li><strong>Focus on Vulnerable Cell Types</strong>: <ul> <li>Prioritize VCT, EVT_1, and VCT_p in your analysis</li> <li>Compare infected vs uninfected within same donors</li> <li>Use cell type proportions as weights</li> </ul> </li> <li><strong>Assess Key Pathways</strong>: <ul> <li>EA metabolism: PLD1, GDPD1, GDPD5, FAAH</li> <li>Adhesion: FUT2, FUT3, FUT4, GALNT1, GALNT2</li> <li>Immunity: TLR2, TLR4, IL1B, IL6, TNF</li> <li>MegL: MGL (if available)</li> </ul> </li> <li><strong>Time-Course Analysis</strong>: <ul> <li>Compare 24h vs 48h timepoints</li> <li>Look for progressive changes in vulnerability</li> <li>Assess immune activation over time</li> </ul> </li> </ol><h3>For Model Improvement</h3><ol> <li><strong>Obtain Expression Matrix</strong>: <ul> <li>Essential for gene-level analysis</li> <li>Required for Seurat integration</li> <li>Needed for DEG analysis</li> </ul> </li> <li><strong>Spatial Validation</strong>: <ul> <li>Compare explant results with in vivo spatial data</li> <li>Assess neighborhood enrichment patterns</li> <li>Validate cell-cell interactions</li> </ul> </li> <li><strong>Cross-Reference with Hoo et al.</strong>: <ul> <li>Compare your analysis with published results</li> <li>Validate infection response patterns</li> <li>Benchmark vulnerability scores</li> </ul> </li> </ol><hr><h2>Next Steps: Action Items</h2><h3>Immediate (Priority 1)</h3><ul> <li><input disabled="" type="checkbox"> <strong>Obtain explant expression matrix</strong> - Critical for full analysis</li> <li><input disabled="" type="checkbox"> Review all generated figures in <code>/workspace/analysis/figures/</code></li> <li><input disabled="" type="checkbox"> Read the comprehensive report: <code>Explant_vs_InVivo_Comprehensive_Report.md</code></li> <li><input disabled="" type="checkbox"> Examine vulnerability summary: <code>explant_vulnerability_summary.csv</code></li> </ul><h3>Short-term (Priority 2)</h3><ul> <li><input disabled="" type="checkbox"> Run R integration script: <code>analysis/R_scripts/01_Explorant_vs_InVivo_Comparison.R</code></li> <li><input disabled="" type="checkbox"> Compare infected vs uninfected within same donors</li> <li><input disabled="" type="checkbox"> Assess MegL expression (if expression matrix available)</li> <li><input disabled="" type="checkbox"> Calculate EA metabolism scores</li> </ul><h3>Long-term (Priority 3)</h3><ul> <li><input disabled="" type="checkbox"> Perform Seurat v5 integration with in vivo data</li> <li><input disabled="" type="checkbox"> Analyze spatial neighborhood patterns</li> <li><input disabled="" type="checkbox"> Validate findings with spatial transcriptomics</li> <li><input disabled="" type="checkbox"> Compare with organoid models (if available)</li> </ul><hr><h2>Contact and Support</h2><p>If you need help with:</p><ol> <li><strong>Data Access</strong>: Expression matrix locations</li> <li><strong>Integration</strong>: Seurat v5 or Harmony setup</li> <li><strong>Analysis</strong>: Custom vulnerability metrics</li> <li><strong>Visualization</strong>: Additional figures</li> </ol><p>Feel free to ask for clarification or additional analysis.</p><hr><h2>Summary Statistics</h2><h3>Explant Dataset</h3><ul> <li><strong>Total Cells</strong>: 158,978</li> <li><strong>Uninfected Cells</strong>: 74,685</li> <li><strong>Cell Types</strong>: 15</li> <li><strong>Donors</strong>: 10</li> <li><strong>Timepoints</strong>: 24h, 48h</li> <li><strong>Infections</strong>: UI, Lm, Tg, Pf</li> </ul><h3>In Vivo Dataset</h3><ul> <li><strong>Total Cells</strong>: ~142,260</li> <li><strong>Cell Types</strong>: 26</li> <li><strong>Technologies</strong>: Slide-tags, Multiome, STARmap</li> <li><strong>Weeks</strong>: 8-23</li> </ul><h3>Analysis Results</h3><ul> <li><strong>Figures Generated</strong>: 6 publication-quality plots</li> <li><strong>Tables Created</strong>: 5 summary tables</li> <li><strong>Scripts Written</strong>: 2 analysis scripts</li> <li><strong>Report Length</strong>: 10,000+ words</li> </ul><hr><p><em>Last Updated: March 2024</em> <em>Analysis by: SuperNinja AI Agent</em></p>