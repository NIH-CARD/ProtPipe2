dat <- data.table::fread("EXAMPLES/basic_example_data/iPSC.csv")
dat_pro <- create_se(dat)

test_that("correctly get PG_counts", {
  pg_counts <- ProtPipe::get_pg_counts(dat_pro)
  expect_s3_class(pg_counts, "data.frame")
  expect_equal(dim(pg_counts), c(42, 2))
})

test_that("correctly plots pg counts", {
  pg_counts <- ProtPipe::plot_pg_counts(dat_pro)
  pg_counts_cond <- ProtPipe::plot_pg_counts(dat_pro, condition = 'base_condition')
  expect_s3_class(pg_counts, "ggplot")
  expect_s3_class(pg_counts_cond, "ggplot")
})

test_that("correctly plots pg intensities", {
  pg_int <- ProtPipe::plot_pg_intensities(dat_pro)
  expect_s3_class(pg_int, "ggplot")
})

test_that("correctly gets CVs", {
  cvs <- ProtPipe::get_CVs(dat_pro, condition = 'base_condition')
  expect_s3_class(cvs, "data.frame")
  expect_equal(dim(cvs), c(57928,3))
})

test_that("correctly plot CVs", {
  p_cvs <- ProtPipe::plot_CVs(dat_pro, condition = 'base_condition')
  expect_s3_class(p_cvs, "ggplot")
})

test_that("correctly gets correlations", {
  cors <- ProtPipe::get_sample_correlation(dat_pro)
  expect_s3_class(cvs, "data.frame")
  expect_equal(dim(cvs), c(1764, 3))
})

test_that("correctly plot sample correlation", {
  p_cors <- ProtPipe::plot_correlation_heatmap(dat_pro)
  expect_s3_class(p_cvs, "ggplot")
})

