# Package index

## Input

- [`protpipe_example_se`](https://nih-card.github.io/ProtPipe2/docs/reference/protpipe_example_se.md)
  : Example SummarizedExperiment for ProtPipe workflows
- [`neuron_differentiation_intensities`](https://nih-card.github.io/ProtPipe2/docs/reference/neuron_differentiation_intensities.md)
  : Bundled intensity table for the neuron differentiation example
- [`neuron_differentiation_metadata`](https://nih-card.github.io/ProtPipe2/docs/reference/neuron_differentiation_metadata.md)
  : Bundled sample metadata for the neuron differentiation example
- [`ipsc_stem_cell_genes`](https://nih-card.github.io/ProtPipe2/docs/reference/ipsc_stem_cell_genes.md)
  : Bundled iPSC stem cell marker genes
- [`create_se()`](https://nih-card.github.io/ProtPipe2/docs/reference/create_se.md)
  : This function takes a data frame of proteomics data and its
  corresponding sample metadata to construct a SummarizedExperiment
  object. It handles detection of intensity columns, validation, and
  synchronization of metadata.
- [`create_se_from_olink()`](https://nih-card.github.io/ProtPipe2/docs/reference/create_se_from_olink.md)
  : Create a SummarizedExperiment Object from Olink NPX Data
- [`create_se_from_soma()`](https://nih-card.github.io/ProtPipe2/docs/reference/create_se_from_soma.md)
  : Create a SummarizedExperiment Object from SomaScan Data
- [`detect_intensity_cols()`](https://nih-card.github.io/ProtPipe2/docs/reference/detect_intensity_cols.md)
  : Identify Numeric Columns by Index
- [`convert_numeric_cols()`](https://nih-card.github.io/ProtPipe2/docs/reference/convert_numeric_cols.md)
  : Safely Convert Character Columns to Numeric Type
- [`olink_all_output()`](https://nih-card.github.io/ProtPipe2/docs/reference/olink_all_output.md)
  : Format Olink NPX Data into Data and Condition Tables
- [`soma_all_output()`](https://nih-card.github.io/ProtPipe2/docs/reference/soma_all_output.md)
  : format soma adat into data and condition dataframes
- [`run_protpipe_shiny()`](https://nih-card.github.io/ProtPipe2/docs/reference/run_protpipe_shiny.md)
  : Run the protpipe Shiny App

## Quality Control

- [`get_pg_counts()`](https://nih-card.github.io/ProtPipe2/docs/reference/get_pg_counts.md)
  : Get Protein Group Counts Per Sample
- [`plot_pg_counts()`](https://nih-card.github.io/ProtPipe2/docs/reference/plot_pg_counts.md)
  : Plot Protein Group Counts
- [`plot_pg_intensities()`](https://nih-card.github.io/ProtPipe2/docs/reference/plot_pg_intensities.md)
  : Plot Boxplots of Sample Intensity Distributions
- [`get_CVs()`](https://nih-card.github.io/ProtPipe2/docs/reference/get_CVs.md)
  : Calculate Coefficient of Variation (CV) for Protein Groups
- [`plot_CVs()`](https://nih-card.github.io/ProtPipe2/docs/reference/plot_CVs.md)
  : Plot Coefficient of Variation (CV) Distributions
- [`get_sample_correlation()`](https://nih-card.github.io/ProtPipe2/docs/reference/get_sample_correlation.md)
  : Calculate Pairwise Sample Correlations
- [`plot_correlation_heatmap()`](https://nih-card.github.io/ProtPipe2/docs/reference/plot_correlation_heatmap.md)
  : Plot a Sample Correlation Heatmap

## Pre-Processing

- [`lod_filter()`](https://nih-card.github.io/ProtPipe2/docs/reference/lod_filter.md)
  : Filter Assay Based on Limit of Detection (LOD)

- [`apply_min_intenisty()`](https://nih-card.github.io/ProtPipe2/docs/reference/apply_min_intenisty.md)
  : Apply Limit of Detection Threshold

- [`filter_unique_proteins()`](https://nih-card.github.io/ProtPipe2/docs/reference/filter_unique_proteins.md)
  : Removes duplicate analytes

- [`filter_overlap()`](https://nih-card.github.io/ProtPipe2/docs/reference/filter_overlap.md)
  :

  Retain proteins present in a specified group of a
  `SummarizedExperiment` object

- [`filter_proteins_by_percent()`](https://nih-card.github.io/ProtPipe2/docs/reference/filter_proteins_by_percent.md)
  : Filter Proteins by Percentage of Valid Values

- [`filter_outlier_samples()`](https://nih-card.github.io/ProtPipe2/docs/reference/filter_outlier_samples.md)
  : Filter Outlier Samples Based on Protein Counts

- [`log2_transform()`](https://nih-card.github.io/ProtPipe2/docs/reference/log2_transform.md)
  : Performs a log2 transform of protein intensity values

- [`z_score()`](https://nih-card.github.io/ProtPipe2/docs/reference/z_score.md)
  : Z-Score Normalization for Proteins Across Samples

- [`mean_normalize()`](https://nih-card.github.io/ProtPipe2/docs/reference/mean_normalize.md)
  : Mean Normalization of Proteomics Data

- [`median_normalize()`](https://nih-card.github.io/ProtPipe2/docs/reference/median_normalize.md)
  : Median Normalization of Proteomics Data

- [`impute()`](https://nih-card.github.io/ProtPipe2/docs/reference/impute.md)
  : Impute Missing Values with a Constant

- [`impute_min()`](https://nih-card.github.io/ProtPipe2/docs/reference/impute_min.md)
  : Impute Missing Values with the Row Minimum

- [`impute_left_dist()`](https://nih-card.github.io/ProtPipe2/docs/reference/impute_left_dist.md)
  : Impute from a Down-Shifted Normal Distribution

- [`batch_correct()`](https://nih-card.github.io/ProtPipe2/docs/reference/batch_correct.md)
  : Correct for Batch Effects

- [`generate_preprocessing_report()`](https://nih-card.github.io/ProtPipe2/docs/reference/generate_preprocessing_report.md)
  : Generate Markdown Report of Preprocessing Steps

- [`has_step()`](https://nih-card.github.io/ProtPipe2/docs/reference/has_step.md)
  : Check if a Processing Step has been Applied

## Clustering / Dimensionality Reduction

- [`get_PCs()`](https://nih-card.github.io/ProtPipe2/docs/reference/get_PCs.md)
  : Calculate Principal Components Analysis (PCA)
- [`plot_PCs()`](https://nih-card.github.io/ProtPipe2/docs/reference/plot_PCs.md)
  : Plot Principal Component Analysis Results
- [`plot_hierarchical_cluster()`](https://nih-card.github.io/ProtPipe2/docs/reference/plot_hierarchical_cluster.md)
  : Plot a Hierarchical Clustering Dendrogram of Samples
- [`get_umap()`](https://nih-card.github.io/ProtPipe2/docs/reference/get_umap.md)
  : Calculate UMAP Dimensionality Reduction
- [`plot_umap()`](https://nih-card.github.io/ProtPipe2/docs/reference/plot_umap.md)
  : Plot UMAP Dimensionality Reduction Results

## Differential Expression

- [`do_limma_binary()`](https://nih-card.github.io/ProtPipe2/docs/reference/do_limma_binary.md)
  : Perform limma differential expression on a SummarizedExperiment
- [`do_t_test_binary()`](https://nih-card.github.io/ProtPipe2/docs/reference/do_t_test_binary.md)
  : Perform t-test differential expression on a SummarizedExperiment
- [`do_anova()`](https://nih-card.github.io/ProtPipe2/docs/reference/do_anova.md)
  : Perform ANOVA-style differential expression across multiple groups
- [`do_comparison_continuous()`](https://nih-card.github.io/ProtPipe2/docs/reference/do_comparison_continuous.md)
  : Perform limma differential expression for a continuous outcome
- [`plot_volcano()`](https://nih-card.github.io/ProtPipe2/docs/reference/plot_volcano.md)
  : Plot a Volcano Plot for Differential Expression Results
- [`plot_correlation_volcano()`](https://nih-card.github.io/ProtPipe2/docs/reference/plot_correlation_volcano.md)
  : Plot a Volcano Plot for Correlation Results
- [`add_entrez()`](https://nih-card.github.io/ProtPipe2/docs/reference/add_entrez.md)
  : Add Entrez Gene IDs to a Data Frame
- [`read_ontology()`](https://nih-card.github.io/ProtPipe2/docs/reference/read_ontology.md)
  : Read a Custom Ontology File
- [`enrich_go()`](https://nih-card.github.io/ProtPipe2/docs/reference/enrich_go.md)
  : Perform Gene Ontology (GO) Over-Representation Analysis (ORA)
- [`enrich_kegg()`](https://nih-card.github.io/ProtPipe2/docs/reference/enrich_kegg.md)
  : Perform KEGG Over-Representation Analysis (ORA)
- [`enrich_terms()`](https://nih-card.github.io/ProtPipe2/docs/reference/enrich_terms.md)
  : Run Over-Representation Analysis for Custom Gene Sets
- [`gse_go()`](https://nih-card.github.io/ProtPipe2/docs/reference/gse_go.md)
  : Perform Gene Ontology (GO) Gene Set Enrichment Analysis (GSEA)
- [`gse_kegg()`](https://nih-card.github.io/ProtPipe2/docs/reference/gse_kegg.md)
  : Perform KEGG Gene Set Enrichment Analysis (GSEA)
- [`gse_terms()`](https://nih-card.github.io/ProtPipe2/docs/reference/gse_terms.md)
  : Run GSEA for Custom Gene Sets
- [`enrich_pathways()`](https://nih-card.github.io/ProtPipe2/docs/reference/enrich_pathways.md)
  : Perform Comprehensive GO and KEGG Pathway Enrichment Analysis

## Abundance Profiling

- [`compare_protein()`](https://nih-card.github.io/ProtPipe2/docs/reference/compare_protein.md)
  : Creates a ggplot bar chart comparing the intensity of a single
  protein either across all samples or grouped by a condition.
- [`plot_proteomics_heatmap()`](https://nih-card.github.io/ProtPipe2/docs/reference/plot_proteomics_heatmap.md)
  : Plot a Proteomics Heatmap
