test_that("limma wrappers run on the bundled example dataset", {
  se <- load_basic_imputed_se()
  treatment_samples <- paste0("Day28_", 1:3)
  control_samples <- paste0("Day0_", 1:3)

  de <- ProtPipe::do_limma(se, treatment_samples = treatment_samples, control_samples = control_samples)
  de_by_condition <- ProtPipe::do_limma_by_condition(
    se,
    condition = "base_condition",
    control_group = "Day0",
    treatment_group = "Day28"
  )

  expect_true(all(c("logFC", "P.Value", "adj.P.Val") %in% names(de)))
  expect_true(all(c("logFC", "P.Value", "adj.P.Val") %in% names(de_by_condition)))
  expect_s3_class(ProtPipe::plot_volcano(de_by_condition, label_col = "PG.Genes"), "ggplot")
})

test_that("binary limma supports covariates and rejects identical groups", {
  meta <- load_basic_metadata(include_batch = TRUE)
  se <- load_basic_imputed_se(meta)

  de <- ProtPipe::do_limma_binary(
    se,
    condition = "base_condition",
    treatment_group = "Day28",
    control_group = "Day0",
    covariates = "batch"
  )

  expect_true(all(c("logFC", "P.Value", "adj.P.Val") %in% names(de)))
  expect_error(
    ProtPipe::do_limma_binary(
      se,
      condition = "base_condition",
      treatment_group = "Day0",
      control_group = "Day0",
      covariates = NULL
    )
  )
})

test_that("continuous comparison works with bundled numeric metadata", {
  meta <- load_basic_metadata(include_numeric = TRUE)
  se <- ProtPipe::impute_min(load_basic_se(meta), 0)

  de <- ProtPipe::do_comparison_continuous(se, "day_num")

  expect_true(all(c("rho", "P.Value", "adj.P.Val") %in% names(de)))
  expect_s3_class(ProtPipe::plot_correlation_volcano(de, label_col = "PG.Genes"), "ggplot")
})

test_that("custom ontology wrappers work on simple local gene sets", {
  universe <- as.character(1:20)
  term2gene <- data.frame(
    term = rep(c("TERM_UP", "TERM_DOWN"), each = 10),
    gene = universe,
    stringsAsFactors = FALSE
  )
  term2name <- data.frame(
    term = c("TERM_UP", "TERM_DOWN"),
    name = c("Up genes", "Down genes"),
    stringsAsFactors = FALSE
  )

  ora <- ProtPipe::enrich_terms(
    gene_id = as.character(1:10),
    all_gene_vector = universe,
    term2gene = term2gene,
    term2name = term2name
  )
  gsea <- ProtPipe::gse_terms(
    gene_list = stats::setNames(c(seq(2, 1.1, length.out = 10), seq(-1.1, -2, length.out = 10)), universe),
    term2gene = term2gene,
    term2name = term2name
  )

  expect_false(is.null(ora))
  expect_false(is.null(gsea))
})

test_that("add_entrez maps gene symbols onto differential results", {
  skip_if_not_installed("org.Hs.eg.db")

  se <- load_basic_imputed_se()
  de <- ProtPipe::do_limma_by_condition(
    se,
    condition = "base_condition",
    control_group = "Day0",
    treatment_group = "Day28"
  )

  mapped <- ProtPipe::add_entrez(de, gene_col = "PG.Genes")

  expect_false(is.null(mapped))
  expect_true("ENTREZID" %in% names(mapped))
  expect_gt(nrow(mapped), 0)
})

test_that("GO and KEGG wrappers handle invalid IDs gracefully", {
  skip_if_not_installed("org.Hs.eg.db")

  invalid_genes <- c("not_a_gene", "still_not_a_gene")
  invalid_ranked <- stats::setNames(c(2, -2), invalid_genes)

  expect_null(ProtPipe::enrich_go(invalid_genes, invalid_genes))
  expect_null(ProtPipe::gse_go(invalid_ranked))
  expect_null(ProtPipe::enrich_kegg(invalid_genes, invalid_genes))
  expect_null(ProtPipe::gse_kegg(invalid_ranked))
})

test_that("enrich_pathways short-circuits on comparisons with no signal", {
  no_signal <- data.frame(
    PG.Genes = c("TP53", "EGFR", "BRCA1"),
    logFC = c(0, 0, 0),
    adj.P.Val = c(1, 1, 1),
    stringsAsFactors = FALSE
  )

  enrichment <- ProtPipe::enrich_pathways(no_signal, gene_col = "PG.Genes")

  expect_true("message" %in% names(enrichment$results))
  expect_match(enrichment$results$message$message[[1]], "no effect-size signal")
})
