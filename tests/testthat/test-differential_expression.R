dat <- data.table::fread("EXAMPLES/basic_example_data/iPSC.csv")
dat_pro <- create_se(dat)

test_that("correctly perform limma comparison", {
  treatment_samples <- c("Day0_1", "Day0_2", "Day0_3")
  control_samples <- c("Day10_1", "Day10_2", "Day10_3")
  dat_pro_imputed <- impute_min(dat_pro, 0)
  dat_pro_imputedr <- impute_left_dist(log2_transform(dat_pro), shift = 5)
  DE <- do_limma(dat_pro_imputed, treatment_samples = treatment_samples, control_samples = control_samples)
  ProtPipe::plot_volcano(DE)

  DE_condition <- do_limma_by_condition(dat_pro_imputed, condition = "base_condition", control_group = "Day0", treatment_group = "Day10")
  ProtPipe::plot_volcano(DE_condition)
})

test_that("correctly perform continuous comparison", {
  dat <- data.table::fread("/Users/jacobepstein/Documents/noberto/tio2_proteomics - Sheet1.csv")
  meta <- data.table::fread("/Users/jacobepstein/Documents/noberto/tio2_metadata - Sheet1 (2).csv")
  dat_pro <- create_se(dat, meta)

  DE <- do_comparison_continuous(dat_pro, "sp.P01027.CO3_MOUSE.Complement.C3.OS.Mus.musculus.OX.10090.GN.C3.PE.1.SV.3")
  ProtPipe::plot_correlation_volcano(DE)

  DE_condition <- do_limma_by_condition(dat_pro_imputed, condition = "base_condition", control_group = "Day0", treatment_group = "Day10")
  ProtPipe::plot_volcano(DE_condition)
})
