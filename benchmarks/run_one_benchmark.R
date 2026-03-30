#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required.")
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript benchmarks/run_one_benchmark.R <benchmark_name>")
}

benchmark_name <- args[[1]]

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
} else if (!requireNamespace("ProtPipe", quietly = TRUE)) {
  stop("Install ProtPipe or pkgload before running this benchmark script.")
}

required_pkgs <- c("SummarizedExperiment", "org.Hs.eg.db", "umap")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(org.Hs.eg.db)
})

example_path <- file.path("EXAMPLES", "basic_example_data", "iPSC.csv")
stem_gene_path <- file.path("EXAMPLES", "basic_example_data", "stem_cell_gene.csv")

raw_dat <- data.table::fread(example_path)
stem_genes <- data.table::fread(stem_gene_path)[[1]]

se_raw <- ProtPipe::create_se(raw_dat)
se_imputed <- ProtPipe::impute_min(se_raw, 0)
de_day0_day28 <- ProtPipe::do_limma_by_condition(
  se_imputed,
  condition = "base_condition",
  control_group = "Day0",
  treatment_group = "Day28"
)

benchmarks <- list(
  se_creation = function() {
    ProtPipe::create_se(raw_dat)
  },
  qc_get_pg_counts = function() {
    ProtPipe::get_pg_counts(se_raw)
  },
  qc_plot_pg_counts = function() {
    ProtPipe::plot_pg_counts(se_raw)
  },
  qc_plot_pg_intensities = function() {
    ProtPipe::plot_pg_intensities(se_raw)
  },
  qc_get_cvs = function() {
    ProtPipe::get_CVs(se_raw, condition = "base_condition")
  },
  qc_plot_cvs = function() {
    ProtPipe::plot_CVs(se_raw, condition = "base_condition")
  },
  qc_get_sample_correlation = function() {
    ProtPipe::get_sample_correlation(se_raw)
  },
  qc_plot_correlation_heatmap = function() {
    ProtPipe::plot_correlation_heatmap(se_raw)
  },
  clustering_get_pcs = function() {
    ProtPipe::get_PCs(se_imputed, condition = "base_condition")
  },
  clustering_plot_pcs = function() {
    ProtPipe::plot_PCs(se_imputed, condition = "base_condition", pc_x = "PC1", pc_y = "PC2")
  },
  clustering_plot_hierarchical = function() {
    ProtPipe::plot_hierarchical_cluster(se_imputed, dist_method = "euclidean", hclust_method = "complete")
  },
  clustering_get_umap = function() {
    set.seed(1)
    ProtPipe::get_umap(se_imputed, condition = "base_condition", neighbors = 5)
  },
  clustering_plot_umap = function() {
    set.seed(1)
    ProtPipe::plot_umap(se_imputed, condition = "base_condition", neighbors = 5)
  },
  dea_day0_vs_day28 = function() {
    ProtPipe::do_limma_by_condition(
      se_imputed,
      condition = "base_condition",
      control_group = "Day0",
      treatment_group = "Day28"
    )
  },
  dea_plot_volcano = function() {
    ProtPipe::plot_volcano(de_day0_day28, label_col = "PG.Genes")
  },
  pathway_go_enrichment = function() {
    ProtPipe::enrich_pathways(
      de_day0_day28,
      gene_col = "PG.Genes",
      go_org = org.Hs.eg.db::org.Hs.eg.db,
      source = "go",
      go_ont = "BP",
      run_kegg = FALSE,
      run_ora = TRUE,
      run_gsea = TRUE
    )
  },
  protein_barchart = function() {
    ProtPipe::compare_protein(
      se_raw,
      "U3KQP1",
      condition = "base_condition",
      selected_groups = c("Day0", "Day28")
    )
  },
  proteomics_heatmap = function() {
    ProtPipe::plot_proteomics_heatmap(
      se_raw,
      protmeta_col = "PG.Genes",
      genes = stem_genes,
      condition = "base_condition",
      cluster_cols = TRUE,
      cluster_rows = TRUE
    )
  },
  total_workflow = function() {
    se <- ProtPipe::create_se(raw_dat)
    ProtPipe::get_pg_counts(se)
    ProtPipe::plot_pg_counts(se)
    ProtPipe::plot_pg_intensities(se)
    ProtPipe::get_CVs(se, condition = "base_condition")
    ProtPipe::plot_CVs(se, condition = "base_condition")
    ProtPipe::get_sample_correlation(se)
    ProtPipe::plot_correlation_heatmap(se)

    se_imp <- ProtPipe::impute_min(se, 0)
    ProtPipe::get_PCs(se_imp, condition = "base_condition")
    ProtPipe::plot_PCs(se_imp, condition = "base_condition", pc_x = "PC1", pc_y = "PC2")
    ProtPipe::plot_hierarchical_cluster(se_imp, dist_method = "euclidean", hclust_method = "complete")
    set.seed(1)
    ProtPipe::get_umap(se_imp, condition = "base_condition", neighbors = 5)
    set.seed(1)
    ProtPipe::plot_umap(se_imp, condition = "base_condition", neighbors = 5)

    de <- ProtPipe::do_limma_by_condition(
      se_imp,
      condition = "base_condition",
      control_group = "Day0",
      treatment_group = "Day28"
    )
    ProtPipe::plot_volcano(de, label_col = "PG.Genes")
    ProtPipe::enrich_pathways(
      de,
      gene_col = "PG.Genes",
      go_org = org.Hs.eg.db::org.Hs.eg.db,
      source = "go",
      go_ont = "BP",
      run_kegg = FALSE,
      run_ora = TRUE,
      run_gsea = TRUE
    )
    ProtPipe::compare_protein(
      se,
      "U3KQP1",
      condition = "base_condition",
      selected_groups = c("Day0", "Day28")
    )
    ProtPipe::plot_proteomics_heatmap(
      se,
      protmeta_col = "PG.Genes",
      genes = stem_genes,
      condition = "base_condition",
      cluster_cols = TRUE,
      cluster_rows = TRUE
    )
  }
)

if (!benchmark_name %in% names(benchmarks)) {
  stop(
    "Unknown benchmark '", benchmark_name, "'. Available benchmarks: ",
    paste(names(benchmarks), collapse = ", ")
  )
}

gc(full = TRUE)
elapsed_sec <- system.time(invisible(benchmarks[[benchmark_name]]()))[["elapsed"]]
cat("benchmark=", benchmark_name, "\n", sep = "")
cat("elapsed_sec=", sprintf("%.6f", elapsed_sec), "\n", sep = "")
