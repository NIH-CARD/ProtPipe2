library(testthat)
setwd("../../..")

dat <- data.table::fread("EXAMPLES/basic_example_data/iPSC.csv")
dat_pro <- create_se(dat)


test_that("correctly filter proteins and samples", {

  dat_pro_no_duplicates <- ProtPipe::filter_unique_proteins(dat_pro, "PG.Genes")
  expect_true(ProtPipe::has_step(dat_pro_filtered, "filter_unique_proteins"))
  expect_equal(nrow(dat_pro) - nrow(dat_pro_no_duplicates), 260)

  dat_pro_filtered <- ProtPipe::filter_proteins_by_percent(dat_pro, 50)
  expect_true(ProtPipe::has_step(dat_pro_filtered, "filter_proteins_by_percent"))
  expect_equal(nrow(dat_pro) - nrow(dat_pro_filtered), 1063)

  dat_pro_ultra_filtered <- ProtPipe::filter_outlier_samples(dat_pro_filtered, sds = 3)
  expect_true(ProtPipe::has_step(dat_pro_ultra_filtered, "filter_outlier_samples"))
  expect_true(ProtPipe::has_step(dat_pro_ultra_filtered, "filter_proteins_by_percent"))
  expect_equal(ncol(dat_pro) - ncol(dat_pro_ultra_filtered), 1)

  dat_common_prots <- ProtPipe::filter_overlap(dat_pro, condition_name = "base_condition")
  expect_true(ProtPipe::has_step(dat_common_prots, "filter_overlap"))
  expect_equal(nrow(dat_pro) - nrow(dat_common_prots), 1458)
})


test_that("test normalization and transformation", {
  dat_z <- ProtPipe::z_score(dat_pro)
  dat_mean_norm <- ProtPipe::mean_normalize(dat_pro)
  dat_median_norm <- ProtPipe::median_normalize(dat_pro)
  log2_dat <- ProtPipe::log2_transform(dat_pro)

  ProtPipe::plot_pg_intensities(dat_pro)
  ProtPipe::plot_pg_intensities(dat_z)
  ProtPipe::plot_pg_intensities(dat_mean_norm)
  ProtPipe::plot_pg_intensities(dat_median_norm)
  ProtPipe::plot_pg_intensities(log2_dat)
})

test_that("test imputation", {
  imputed_fixed <- ProtPipe::impute(dat_pro, 5)
  freq_table <- table(assay(imputed_fixed))
  mode <- as.numeric(names(freq_table)[which.max(freq_table)])
  expect_equal(mode, 5)
  expect_true(ProtPipe::has_step(imputed_fixed, "impute"))
  imputed_min <- ProtPipe::impute_min(dat_pro, 0.5)
  expect_true(ProtPipe::has_step(imputed_min, "impute_min"))
  imputed_dist <- ProtPipe::impute_left_dist(dat_pro)
  expect_true(ProtPipe::has_step(imputed_dist, "impute_left_dist"))
})

test_that("test batch correction", {
  vdat <- data.table::fread("EXAMPLES/VIRUS/virus_data.tsv")
  vmeta <- data.table::fread("EXAMPLES/VIRUS/virus_metadata.tsv")
  dat_v <- create_se(vdat, sample_metadata = vmeta) %>% impute(0)
  plot_PCs(dat_v, condition = "viral.exposure")
  dat_v_corrected <- ProtPipe::batch_correct(dat_v, batch_variable = "viral.exposure")
  plot_PCs(dat_v_corrected, condition = "viral.exposure")
  dat_v_corrected_careful <- ProtPipe::batch_correct(dat_v, batch_variable = "viral.exposure",
                                                     bio_variables = c("concentration", "time"))
  plot_PCs(dat_v_corrected_careful, condition = "viral.exposure")
})

