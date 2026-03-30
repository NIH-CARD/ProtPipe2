test_that("protein group count helpers return expected shapes", {
  se <- load_basic_se()

  pg_counts <- ProtPipe::get_pg_counts(se)
  expect_s3_class(pg_counts, "data.frame")
  expect_equal(dim(pg_counts), c(ncol(se), 2))

  expect_s3_class(ProtPipe::plot_pg_counts(se), "ggplot")
  expect_s3_class(ProtPipe::plot_pg_counts(se, condition = "base_condition"), "ggplot")
})

test_that("intensity and CV QC plots run on the example dataset", {
  se <- load_basic_se()

  cvs <- ProtPipe::get_CVs(se, condition = "base_condition")
  expect_s3_class(cvs, "data.frame")
  expect_true(all(c("Protein", "CV", "Condition") %in% names(cvs)))

  expect_s3_class(ProtPipe::plot_pg_intensities(se), "ggplot")
  expect_s3_class(ProtPipe::plot_CVs(se, condition = "base_condition"), "ggplot")
})

test_that("sample correlation helpers return expected shapes", {
  se <- load_basic_se()

  cors <- ProtPipe::get_sample_correlation(se)
  expect_s3_class(cors, "data.frame")
  expect_equal(dim(cors), c(ncol(se) ^ 2, 3))

  p <- ProtPipe::plot_correlation_heatmap(se)
  expect_s3_class(p, "ggplot")
})
