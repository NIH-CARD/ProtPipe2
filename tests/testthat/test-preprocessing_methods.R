test_that("threshold-based filtering records processing steps", {
  se <- load_basic_se()
  filtered <- ProtPipe::apply_min_intenisty(se, 1000)

  expect_true(ProtPipe::has_step(filtered, "apply_min_intenisty"))
  expect_gt(sum(is.na(SummarizedExperiment::assay(filtered))),
            sum(is.na(SummarizedExperiment::assay(se))))
})

test_that("LOD filtering masks values below the per-protein threshold", {
  se <- load_basic_se()
  SummarizedExperiment::rowData(se)$Buffer <- matrixStats::rowMeans2(
    SummarizedExperiment::assay(se),
    na.rm = TRUE
  )

  filtered <- ProtPipe::lod_filter(se, lod_col = "Buffer")

  expect_true(ProtPipe::has_step(filtered, "lod_filter"))
  expect_gt(sum(is.na(SummarizedExperiment::assay(filtered))),
            sum(is.na(SummarizedExperiment::assay(se))))
})

test_that("row and sample filters preserve object validity", {
  se <- load_basic_se()

  deduped <- ProtPipe::filter_unique_proteins(se, "PG.Genes")
  percent_filtered <- ProtPipe::filter_proteins_by_percent(se, 50)
  outlier_filtered <- ProtPipe::filter_outlier_samples(percent_filtered, sds = 3)
  overlap_filtered <- ProtPipe::filter_overlap(se, condition_name = "base_condition")

  expect_true(ProtPipe::has_step(deduped, "filter_unique_proteins"))
  expect_true(ProtPipe::has_step(percent_filtered, "filter_proteins_by_percent"))
  expect_true(ProtPipe::has_step(outlier_filtered, "filter_outlier_samples"))
  expect_true(ProtPipe::has_step(overlap_filtered, "filter_overlap"))

  expect_lte(nrow(deduped), nrow(se))
  expect_lte(nrow(percent_filtered), nrow(se))
  expect_lte(nrow(overlap_filtered), nrow(se))
  expect_lte(ncol(outlier_filtered), ncol(percent_filtered))
})

test_that("normalization and transformation methods log their processing steps", {
  se <- load_basic_se()

  z <- ProtPipe::z_score(se)
  mean_norm <- ProtPipe::mean_normalize(se)
  median_norm <- ProtPipe::median_normalize(se)
  log2_se <- ProtPipe::log2_transform(se)

  expect_true(ProtPipe::has_step(z, "z_score"))
  expect_true(ProtPipe::has_step(mean_norm, "mean_normalize"))
  expect_true(ProtPipe::has_step(median_norm, "median_normalize"))
  expect_true(ProtPipe::has_step(log2_se, "log2_transform"))
  expect_equal(dim(z), dim(se))
  expect_equal(dim(mean_norm), dim(se))
  expect_equal(dim(median_norm), dim(se))
  expect_equal(dim(log2_se), dim(se))
})

test_that("imputation methods record their processing steps", {
  se <- load_basic_se()

  imputed_fixed <- ProtPipe::impute(se, 5)
  imputed_min <- ProtPipe::impute_min(se, 0.5)
  imputed_dist <- ProtPipe::impute_left_dist(se)

  expect_true(ProtPipe::has_step(imputed_fixed, "impute"))
  expect_true(ProtPipe::has_step(imputed_min, "impute_min"))
  expect_true(ProtPipe::has_step(imputed_dist, "impute_left_dist"))
  expect_false(anyNA(SummarizedExperiment::assay(imputed_fixed)))
})

test_that("batch correction works with bundled example metadata", {
  meta <- load_basic_metadata(include_batch = TRUE)
  se <- load_basic_se(meta)
  corrected <- ProtPipe::batch_correct(ProtPipe::impute(se, 0), batch_variable = "batch")

  expect_true(ProtPipe::has_step(corrected, "batch_corrected"))
  expect_equal(dim(SummarizedExperiment::assay(corrected)),
               dim(SummarizedExperiment::assay(se)))
})

test_that("generate_preprocessing_report writes a markdown file", {
  se <- ProtPipe::log2_transform(load_basic_se())
  report <- tempfile(fileext = ".md")

  out <- ProtPipe::generate_preprocessing_report(se, output_file = report)

  expect_equal(out, report)
  expect_true(file.exists(report))
  expect_match(paste(readLines(report), collapse = "\n"), "log2_transform")
})
