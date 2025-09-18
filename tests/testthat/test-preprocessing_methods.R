library(testthat)
setwd("../../..")

dat <- data.table::fread("EXAMPLES/basic_example_data/iPSC.csv")
dat_pro <- create_protdata(dat)


test_that("correctly filter proteins and samples", {
  dat_pro_filtered <- ProtPipe::filter_proteins_by_percent(dat_pro, 50)
  expect_true(ProtPipe::has_step(dat_pro_filtered, "filter_proteins_by_percent"))
  dat_pro_ultra_filtered <- ProtPipe::filter_outlier_samples(dat_pro_filtered, sds = 2)
  expect_true(ProtPipe::has_step(dat_pro_ultra_filtered, "filter_outlier_samples"))
})


test_that("test full pipeline", {
  df <- data.table::fread("EXAMPLES/DIANN/report.pg_matrix.tsv")
  dat_pro <- create_protdata(df, intensity_cols = c(5:138))
  expect_s4_class(dat_pro, "ProtData")
})

