# Lucas Lab BLA snRNA-seq Analysis Codebase

This repository contains the analysis scripts associated with the Teravskis, Baumgartner, and Lucas 2026 basolateral amygdala (BLA) single-nucleus RNA-seq study. The code supports preprocessing, quality control, clustering, cell-type annotation, barcode reassignment, differential expression, enrichment analysis, regulatory-network analysis, spatial/Tangram analysis, and figure generation.

This is primarily a script-only code release. Raw data, intermediate data objects, live-data test inputs, generated figures, and local test logs are not included in this upload folder. The only bundled dependency artifacts are a local Tangram fork and a local scvi-tools drop in patch used by the spatial-mapping and reintegration/reassignment workflows respectively. The scripts are intended to document and support reproducible execution of the analytical workflow used for the paper, subject to availability of the required input data and the appropriate R/Python package environment.

## Repository Structure

```text
analyses/
  cell_type_annotation/     Summary tables for MapMyCells annotations.
  comparative/              RRHO2 and cross-condition DEG comparisons.
  differential_expression/  Seurat marker and cluster-wise DEG workflows.
  enrichment/               GO GSEA/ORA workflows.
  qc/                       Barcode, sample, and integration-quality summaries.
  regulatory/               DoRothEA/decoupleR TF activity and regulon networks.

pipelines/
  cell_type_annotation/     Garnett classifier workflow.
  clustering/               SCTransform clustering passes.
  integration/              Seurat reintegration template.
  network_linearization/    STRING/GO network linearization workflow.
  preprocessing/            10x HDF5 filtering, SoupX, and doublet removal.
  reassignment/             H5AD export, scVI/scANVI, barcode reassignment, and related tests.
  spatial_mapping/tangram/  Tangram alignment, spatial summaries, and slice plots.

utilities/
  General conversion, barcode counting, HDF5 merging, and DEG annotation helpers.

visualization/
  cluster/                  Cluster-colored UMAP plots.
  deg/                      DEG volcano, dot, Venn, and heatmap plots.
  enrichment/               GSEA/ORA count, heatmap, and top-term plots.
  expression/               Feature, receptor, rug/ridge, and dendrogram plots.
  general/                  General-purpose plotting helpers.
  mapmycells/               MapMyCells stacked bars, pie charts, UMAPs, and doublet plots.
  regulatory/               TF activity and top-TF summary plots.

forks/
  tangram_modified/         Local installable Tangram v1.0.4 fork used by Tangram alignment.
```

## Command-line Refactor

The scripts in this release were refactored from an iterative local analysis workflow into command-line driven scripts. Most R scripts can be inspected with:

```bash
Rscript path/to/script.R --help
```

Python scripts can generally be inspected with:

```bash
python path/to/script.py --help
```

The refactor standardizes common command-line patterns: explicit input paths, explicit output paths, configurable thresholds, reproducible seeds where appropriate, and output directories that can be regenerated without editing hard-coded local paths. Some scripts still represent workflow templates or specialized paper-analysis utilities and may require input objects with study-specific columns or metadata.

## Tangram Fork and License

The Tangram alignment workflow uses a local forked package staged under:

```text
forks/tangram_modified/
```

The installable source distribution is:

```text
forks/tangram_modified/dist/tangram_modified-1.0.4+pjtfork.tar.gz
```

This fork is based on Tangram v1.0.4 from the Broad Institute Tangram repository and was modified to expose diagnostic training-loss curves, including total loss, main loss, and KL regularization loss. It is used by `pipelines/spatial_mapping/tangram/tangram_alignment.py` through the standard `import tangram as tg` import path after installation.

The original Tangram project is distributed under the BSD 3-Clause License. The visible license file for the fork is included at:

```text
forks/tangram_modified/LICENSE.lic
```

That file preserves the Broad Institute copyright notice, the BSD 3-Clause license text, and the modification notice for this local fork. The original Tangram source repository is: https://github.com/broadinstitute/Tangram

To install the fork from this staged repository:

```bash
pip install forks/tangram_modified/dist/tangram_modified-1.0.4+pjtfork.tar.gz
```

## scvi-tools Drop In Patch and License

The scvi-tools workflow uses a local drop in patch of scvi-tools v1.3.3 extracted under:

```text
forks/scvi_tools_drop_in
```

The associated file perserves the Yosef Lab, Weizmann Institute of Science copyright notice, the BSD 3-Clause license text, and the modification notice for this local patch. The original scvi-tools source repository is: https://github.com/scverse/scvi-tools

The installation instructions are found at:

```text
forks/scvi_tools_drop_in/INSTALL.md
```

## Testing Status

During refactor, the codebase was tested using a combination of:

- CLI smoke tests to confirm scripts parse arguments, load dependencies, and expose help text.
- Minimal fixture tests to confirm scripts can run on small shaped inputs.
- Live-data tests on available study data to confirm outputs are generated in the expected formats.
- Output verification for selected files, including checking CSV/table structure, PDF/PNG generation, non-empty output files, and visual rendering checks for plots where layout mattered.

The test process was intended to catch broken paths, missing arguments, dependency failures, parser mismatches, empty outputs, and plot-generation failures introduced during CLI refactoring. Several scripts were updated during this process to make command-line execution more robust and to better support downstream script chaining.

## Testing Limitations

Testing cannot guarantee that every possible parameter combination, package version, operating system, or input schema variation has been exercised. Some workflows are computationally intensive and depend on large Seurat, AnnData, 10x HDF5, spatial transcriptomics, or enrichment-result objects that are not included in this script-only release. Full reproduction also depends on external packages, online resources, reference databases, and data objects that may change over time.

Users should treat this repository as paper-associated analysis code rather than a general-purpose software package. Before reusing the scripts on new data, verify that input columns, identifiers, factor levels, genome/species assumptions, and output paths match the intended workflow. For publication or reuse, rerun scripts in a controlled environment and inspect generated tables and figures directly.

## Script Inventory

### Analyses

| Script | Purpose |
| --- | --- |
| `analyses/cell_type_annotation/summarize_mapmycells_annotations.R` | Summarise MapMyCells hierarchical annotations by UMAP/Seurat cluster. |
| `analyses/comparative/rrho2_cluster_overlap_heatmaps.R` | Generate RRHO2 heatmaps comparing ranked cluster-wise DEG signatures across condition contrasts. |
| `analyses/differential_expression/run_clusterwise_differential_expression.R` | Run cluster-wise Seurat differential expression across specified group comparisons. |
| `analyses/differential_expression/run_findallmarkers_on_seurat.R` | Run Seurat `FindAllMarkers` on a serialized Seurat object and export the marker table. |
| `analyses/differential_expression/summarize_deg_overlap_counts.R` | Summarise overlapping positive and negative DEGs across pairwise comparisons. |
| `analyses/enrichment/cluster_profiler/by_cluster_gsea_wang_similarity.R` | Run cluster-wise GO GSEA with clusterProfiler and export ontology-level results. |
| `analyses/enrichment/cluster_profiler/by_cluster_ora_wang_similarity.R` | Run cluster-wise GO ORA with optional Wang semantic-similarity simplification. |
| `analyses/enrichment/run_cluster_gsea_core_enrichment.R` | Run cluster-wise GO GSEA and export core-enrichment gene information. |
| `analyses/qc/compute_lisi_scores.R` | Compute LISI integration metrics from Seurat embeddings and export diagnostic plots. |
| `analyses/qc/summarize_barcodes_by_cluster_and_sample.R` | Count barcodes per cluster and sample prefix. |
| `analyses/regulatory/regulon_network_omnibus.R` | Build DoRothEA/decoupleR-derived TF-target networks with optional STRING expansion. |
| `analyses/regulatory/score_tf_activity_dorothea_decouplr.R` | Score TF activity from DEG tables using DoRothEA regulons and decoupleR WMEAN or MLM methods. |

### Pipelines

| Script | Purpose |
| --- | --- |
| `pipelines/cell_type_annotation/run_garnett_cluster_classifier.R` | Train and apply a Garnett classifier using selected training and target clusters. |
| `pipelines/clustering/merged_subsets_sctransform.R` | Merge multiple 10x HDF5 subsets, run SCTransform, and export clustering outputs. |
| `pipelines/clustering/step01_first_pass_sctransform.R` | Run first-pass SCTransform and variable-feature diagnostics from a 10x HDF5 matrix. |
| `pipelines/clustering/step02_second_pass_sctransform.R` | Run second-pass SCTransform, PCA, UMAP, clustering, and marker export using selected features. |
| `pipelines/integration/seurat_cluster_integration_template.R` | Template for reintegrating a cluster of interest against the remaining Seurat object. |
| `pipelines/network_linearization/make_network_graph_colorized.R` | Render a color-coded STRING network using GO molecular-function annotations. |
| `pipelines/network_linearization/step01_unique_genes_pfilter_multicsv.R` | Collect unique significant DEG genes across multiple CSV files. |
| `pipelines/network_linearization/step02_string_network_plot.R` | Retrieve STRING interactions and export a network plot and adjacency matrix. |
| `pipelines/network_linearization/step03_cmds_sort_by_connection.R` | Use classical multidimensional scaling to order genes from a similarity matrix. |
| `pipelines/network_linearization/step04_convert_string_id_to_gene_symbol.R` | Convert STRING IDs to preferred gene symbols. |
| `pipelines/network_linearization/step05_basic_ora_wang_similarity.R` | Run GO ORA on network genes with optional Wang semantic simplification. |
| `pipelines/network_linearization/step06_create_ranked_heatmap.R` | Create DEG heatmaps ordered by STRING-derived gene ranks. |
| `pipelines/preprocessing/filter_high_mt_10x_h5.R` | Remove high-mitochondrial barcodes from a 10x HDF5 matrix. |
| `pipelines/preprocessing/filter_low_umi_clusters_10x_h5.R` | Remove specified low-quality clusters from a 10x HDF5 matrix. |
| `pipelines/preprocessing/remove_doublets_scdblfinder_10x_h5.R` | Remove doublet-classified barcodes from a 10x HDF5 matrix. |
| `pipelines/preprocessing/run_scdblfinder_on_seurat.R` | Batch-run scDblFinder on Seurat RDS files and export doublet/QC summaries. |
| `pipelines/preprocessing/soupx_manual_correction.R` | Run SoupX background correction using filtered/raw 10x matrices and optional annotations. |
| `pipelines/reassignment/chi_squared_from_csv.R` | Compute and visualize chi-squared residuals from CSV/count-table inputs. |
| `pipelines/reassignment/create_h5ad.R` | Split a Seurat object into training/target subsets and export H5AD files. |
| `pipelines/reassignment/differential_expression_complete.R` | Run differential-expression summaries for barcode reassignment workflows. |
| `pipelines/reassignment/make_umap_after_reintegration.R` | Generate UMAP visualizations after reintegration and barcode reassignment. |
| `pipelines/reassignment/reassign_barcodes.R` | Reassign selected barcodes using prediction tables and cluster metadata. |
| `pipelines/reassignment/run_create_scvimodel_kwargs.py` | Train or reuse SCVI/SCANVI models and export reassignment predictions. |
| `pipelines/reassignment/stacked_bar_chi_squared_from_count_table.R` | Build stacked bar and chi-squared summaries from count tables. |
| `pipelines/spatial_mapping/tangram/allen_region_barplot.R` | Summarise Tangram assignments across Allen Brain Atlas parcellation substructures. |
| `pipelines/spatial_mapping/tangram/custom_hex_merfish_slice_plots.R` | Generate per-slice MERFISH scatter plots using custom cluster colors. |
| `pipelines/spatial_mapping/tangram/sliceplot_cluster_probability_coronal.py` | Plot coronal slice-level Tangram cluster-probability maps. |
| `pipelines/spatial_mapping/tangram/sliceplot_cluster_probability_sagittal.py` | Plot sagittal slice-level Tangram cluster-probability maps. |
| `pipelines/spatial_mapping/tangram/step01_bin_by_ccf_coordinates.R` | Bin Tangram spatial coordinates and aggregate cell counts by cluster. |
| `pipelines/spatial_mapping/tangram/step02_coordinate_line_plots_loess.R` | Plot binned spatial cluster distributions with optional LOESS smoothing. |
| `pipelines/spatial_mapping/tangram/tangram_alignment.py` | Run Tangram alignment between single-cell and spatial AnnData objects. |

### Utilities

| Script | Purpose |
| --- | --- |
| `utilities/annotate_deg.R` | Add NCBI/organism database gene descriptions to DEG tables. |
| `utilities/convert_10x_h5_to_h5ad.R` | Convert 10x HDF5 matrices to AnnData H5AD files. |
| `utilities/count_barcodes_from_10x_h5.R` | Count barcodes in a 10x HDF5 matrix. |
| `utilities/count_barcodes_per_cluster.R` | Count barcodes per cluster and optional sample prefix from CSV metadata. |
| `utilities/merge_10x_h5_matrices.R` | Merge multiple 10x HDF5 matrices into a combined HDF5 matrix. |

### Visualization

| Script | Purpose |
| --- | --- |
| `visualization/cluster/plot_cluster_colored_umap.R` | Plot UMAP coordinates colored by a cluster-to-color lookup table. |
| `visualization/deg/deg_volcano_plot.R` | Generate DEG volcano plots from tabular differential-expression results. |
| `visualization/deg/plot_deg_comparison_venn_diagrams.R` | Create scaled DEG Venn diagrams across pairwise comparisons and clusters. |
| `visualization/deg/plot_deg_counts_stacked_bar.R` | Plot stacked DEG counts by cluster, comparison, and direction. |
| `visualization/deg/plot_deg_dot_manual_selection.R` | Generate DEG dot plots for manually selected genes. |
| `visualization/deg/plot_deg_volcano_manual_selection.R` | Generate volcano plots for manually selected DEG sets or genes. |
| `visualization/deg/plot_findallmarkers_dot_manual.R` | Generate dot plots from `FindAllMarkers` output using manual gene lists. |
| `visualization/deg/plot_findallmarkers_volcano.R` | Generate volcano plots for all clusters in a `FindAllMarkers` table. |
| `visualization/deg/plot_findallmarkers_volcano_manual.R` | Generate cluster/gene-specific volcano plots from `FindAllMarkers` output. |
| `visualization/deg/plot_log2fc_complex_heatmap.R` | Create ComplexHeatmap visualizations from log2FC matrices. |
| `visualization/enrichment/plot_gsea_cluster_heatmap.R` | Generate wide-format GSEA summary tables and clustered heatmaps. |
| `visualization/enrichment/plot_gsea_term_count_histograms.R` | Plot enriched/de-enriched GSEA term counts across clusters. |
| `visualization/enrichment/plot_ora_term_count_histograms.R` | Plot ORA term counts across clusters and comparisons. |
| `visualization/enrichment/plot_top_go_terms_upset.R` | Produce top GO-term dot/bar summaries and optional UpSet diagrams. |
| `visualization/expression/loglog_hormone_receptor_plots.R` | Generate hormone-receptor pseudobulk log-log plots and DEG summaries. |
| `visualization/expression/plot_expression_rug_and_ridge.R` | Generate rug and ridge plots for ORA/GSEA gene-set positions. |
| `visualization/expression/plot_feature_dendrogram.R` | Compute and plot cluster dendrograms from feature expression. |
| `visualization/expression/plot_feature_expression.R` | Generate Seurat feature plots for selected genes. |
| `visualization/expression/plot_hormone_receptor_bubbles.R` | Generate hormone-receptor bubble plots from DEG tables. |
| `visualization/expression/plot_variable_features_heatmap.R` | Plot selected variable-feature expression using ComplexHeatmap. |
| `visualization/general/plot_dot_from_csv.R` | General-purpose dot plot from a CSV table. |
| `visualization/mapmycells/mapmycells_doublet_comparison_plot.R` | Compare MapMyCells doublet estimates with DoubletFinder/scDblFinder scores. |
| `visualization/mapmycells/plot_mapmycells_pies_and_umap.R` | Generate MapMyCells pie charts and annotation-colored UMAPs. |
| `visualization/mapmycells/plot_mapmycells_stacked_bars.R` | Generate percentage stacked bars for MapMyCells classifications. |
| `visualization/regulatory/plot_tf_activity_bars.R` | Plot TF activity scores across comparisons and clusters. |
| `visualization/regulatory/plot_top_tfs_upset.R` | Summarise top TFs across clusters with dot/bar and optional UpSet plots. |

## General Usage

A typical run pattern is:

```bash
Rscript analyses/regulatory/score_tf_activity_dorothea_decouplr.R \
  --deg-csv path/to/clusterwise_deg.csv \
  --output-dir path/to/output/tf_activity \
  --method wmean \
  --output-format both
```

or:

```bash
python pipelines/spatial_mapping/tangram/tangram_alignment.py \
  --sc-h5ad path/to/single_cell_reference.h5ad \
  --sp-h5ad path/to/spatial_data.h5ad \
  --gene-map-csv path/to/gene_map.csv \
  --output-dir path/to/output/tangram
```

Use each script's `--help` output for the most specific input schema and options.

## Citation

If using this code, please cite the associated Teravskis, Baugartner, and Lucas 2026 paper.
