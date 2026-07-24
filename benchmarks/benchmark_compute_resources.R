#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required.")
  if (!requireNamespace("peakRAM", quietly = TRUE)) stop("Package 'peakRAM' is required.")
})

if (!requireNamespace("ProtPipe", quietly = TRUE)) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Install ProtPipe or pkgload before running this benchmark script.")
  }
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
}

required_pkgs <- c("SummarizedExperiment", "org.Hs.eg.db")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

iterations <- 3L
out_dir <- file.path("benchmarks", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

example_path <- file.path("EXAMPLES", "basic_example_data", "iPSC.csv")
stem_gene_path <- file.path("EXAMPLES", "basic_example_data", "stem_cell_gene.csv")

if (!file.exists(example_path)) stop("Missing example dataset: ", example_path)
if (!file.exists(stem_gene_path)) stop("Missing heatmap gene list: ", stem_gene_path)

extract_peakram_metrics <- function(x) {
  elapsed_col <- intersect(c("Elapsed_Time_sec", "Elapsed_Time_seconds"), names(x))[1]
  peak_col <- intersect(c("Peak_RAM_Used_MiB", "Peak_RAM_Used_MB"), names(x))[1]

  if (is.na(elapsed_col) || is.na(peak_col)) {
    stop("Unexpected peakRAM output columns: ", paste(names(x), collapse = ", "))
  }

  list(
    elapsed_sec = as.numeric(x[[elapsed_col]][1]),
    peak_ram_mib = as.numeric(x[[peak_col]][1])
  )
}

bench_expr <- function(name, expr, iterations = 3L, env = parent.frame()) {
  expr <- substitute(expr)

  runs <- lapply(seq_len(iterations), function(i) {
    gc(full = TRUE)
    res <- peakRAM::peakRAM(eval(expr, env))
    metrics <- extract_peakram_metrics(res)
    data.frame(
      benchmark = name,
      iteration = i,
      runtime_sec = metrics$elapsed_sec,
      peak_ram_mib = metrics$peak_ram_mib,
      stringsAsFactors = FALSE
    )
  })

  raw <- do.call(rbind, runs)
  summary <- data.frame(
    benchmark = name,
    iterations = iterations,
    median_runtime_sec = stats::median(raw$runtime_sec),
    peak_ram_mib = max(raw$peak_ram_mib),
    stringsAsFactors = FALSE
  )

  list(raw = raw, summary = summary)
}

raw_dat <- data.table::fread(example_path)
stem_genes <- data.table::fread(stem_gene_path)[[1]]

se_raw <- ProtPipe2::create_se(raw_dat)
se_imputed <- ProtPipe2::impute_min(se_raw, 0)
de_day0_day28 <- ProtPipe2::do_limma_by_condition(
  se_imputed,
  condition = "base_condition",
  control_group = "Day0",
  treatment_group = "Day28"
)

benchmarks <- list(
  bench_expr("se_creation", ProtPipe2::create_se(raw_dat), iterations),
  bench_expr("qc_get_pg_counts", ProtPipe2::get_pg_counts(se_raw), iterations),
  bench_expr("qc_plot_pg_counts", ProtPipe2::plot_pg_counts(se_raw), iterations),
  bench_expr("qc_plot_pg_intensities", ProtPipe2::plot_pg_intensities(se_raw), iterations),
  bench_expr("qc_get_cvs", ProtPipe2::get_CVs(se_raw, condition = "base_condition"), iterations),
  bench_expr("qc_plot_cvs", ProtPipe2::plot_CVs(se_raw, condition = "base_condition"), iterations),
  bench_expr("qc_get_sample_correlation", ProtPipe2::get_sample_correlation(se_raw), iterations),
  bench_expr("qc_plot_correlation_heatmap", ProtPipe2::plot_correlation_heatmap(se_raw), iterations),
  bench_expr("clustering_get_pcs", ProtPipe2::get_PCs(se_imputed, condition = "base_condition"), iterations),
  bench_expr(
    "clustering_plot_pcs",
    ProtPipe2::plot_PCs(se_imputed, condition = "base_condition", pc_x = "PC1", pc_y = "PC2"),
    iterations
  ),
  bench_expr(
    "clustering_plot_hierarchical",
    ProtPipe2::plot_hierarchical_cluster(se_imputed, dist_method = "euclidean", hclust_method = "complete"),
    iterations
  ),
  bench_expr(
    "clustering_get_umap",
    {
      set.seed(1)
      ProtPipe2::get_umap(se_imputed, condition = "base_condition", neighbors = 5)
    },
    iterations
  ),
  bench_expr(
    "clustering_plot_umap",
    {
      set.seed(1)
      ProtPipe2::plot_umap(se_imputed, condition = "base_condition", neighbors = 5)
    },
    iterations
  ),
  bench_expr(
    "dea_day0_vs_day28",
    ProtPipe2::do_limma_by_condition(
      se_imputed,
      condition = "base_condition",
      control_group = "Day0",
      treatment_group = "Day28"
    ),
    iterations
  ),
  bench_expr(
    "dea_plot_volcano",
    ProtPipe2::plot_volcano(de_day0_day28, label_col = "PG.Genes"),
    iterations
  ),
  bench_expr(
    "pathway_go_enrichment",
    ProtPipe2::enrich_pathways(
      de_day0_day28,
      gene_col = "PG.Genes",
      go_org = org.Hs.eg.db::org.Hs.eg.db,
      source = "go",
      go_ont = "BP",
      run_kegg = FALSE,
      run_ora = TRUE,
      run_gsea = TRUE
    ),
    iterations
  ),
  bench_expr(
    "protein_barchart",
    ProtPipe2::compare_protein(
      se_raw,
      "U3KQP1",
      condition = "base_condition",
      selected_groups = c("Day0", "Day28")
    ),
    iterations
  ),
  bench_expr(
    "proteomics_heatmap",
    ProtPipe2::plot_proteomics_heatmap(
      se_raw,
      protmeta_col = "PG.Genes",
      genes = stem_genes,
      condition = "base_condition",
      cluster_cols = TRUE,
      cluster_rows = TRUE
    ),
    iterations
  ),
  bench_expr(
    "total_workflow",
    {
      se <- ProtPipe2::create_se(raw_dat)
      ProtPipe2::get_pg_counts(se)
      ProtPipe2::plot_pg_counts(se)
      ProtPipe2::plot_pg_intensities(se)
      ProtPipe2::get_CVs(se, condition = "base_condition")
      ProtPipe2::plot_CVs(se, condition = "base_condition")
      ProtPipe2::get_sample_correlation(se)
      ProtPipe2::plot_correlation_heatmap(se)

      se_imp <- ProtPipe2::impute_min(se, 0)
      ProtPipe2::get_PCs(se_imp, condition = "base_condition")
      ProtPipe2::plot_PCs(se_imp, condition = "base_condition", pc_x = "PC1", pc_y = "PC2")
      ProtPipe2::plot_hierarchical_cluster(se_imp, dist_method = "euclidean", hclust_method = "complete")
      set.seed(1)
      ProtPipe2::get_umap(se_imp, condition = "base_condition", neighbors = 5)
      set.seed(1)
      ProtPipe2::plot_umap(se_imp, condition = "base_condition", neighbors = 5)

      de <- ProtPipe2::do_limma_by_condition(
        se_imp,
        condition = "base_condition",
        control_group = "Day0",
        treatment_group = "Day28"
      )
      ProtPipe2::plot_volcano(de, label_col = "PG.Genes")
      ProtPipe2::enrich_pathways(
        de,
        gene_col = "PG.Genes",
        go_org = org.Hs.eg.db::org.Hs.eg.db,
        source = "go",
        go_ont = "BP",
        run_kegg = FALSE,
        run_ora = TRUE,
        run_gsea = TRUE
      )
      ProtPipe2::compare_protein(
        se,
        "U3KQP1",
        condition = "base_condition",
        selected_groups = c("Day0", "Day28")
      )
      ProtPipe2::plot_proteomics_heatmap(
        se,
        protmeta_col = "PG.Genes",
        genes = stem_genes,
        condition = "base_condition",
        cluster_cols = TRUE,
        cluster_rows = TRUE
      )
    },
    iterations
  )
)

summary_df <- do.call(rbind, lapply(benchmarks, `[[`, "summary"))
raw_df <- do.call(rbind, lapply(benchmarks, `[[`, "raw"))

summary_path <- file.path(out_dir, paste0("compute_benchmark_summary_", timestamp, ".csv"))
raw_path <- file.path(out_dir, paste0("compute_benchmark_raw_", timestamp, ".csv"))

utils::write.csv(summary_df, summary_path, row.names = FALSE)
utils::write.csv(raw_df, raw_path, row.names = FALSE)

cat("Dataset:", example_path, "\n")
cat("Proteins:", nrow(se_raw), "Samples:", ncol(se_raw), "\n")
cat("Iterations per benchmark:", iterations, "\n")
cat("Summary results:", summary_path, "\n")
cat("Raw results:", raw_path, "\n\n")
print(summary_df, row.names = FALSE)
