dat <- SomaDataIO::read_adat(package_path("EXAMPLES", "soma", "example_data_v5.0_plasma.adat"))
data1 <- soma_all_output(dat)


test_that("create se object from somascan data without filtering", {
  soma_pro <- create_se_from_soma(dat, filter = FALSE)
  soma_imputed <- impute(soma_pro, 0)
  DE <- do_limma_by_condition(soma_imputed, condition = "Sex", control_group = "F", treatment_group = "M")
  p <- plot_volcano(DE, label_col = "Genes")
  expect_s4_class(soma_pro, "SummarizedExperiment")
})

test_that("create se object from somascan data with filtering", {
  soma_pro <- create_se_from_soma(dat, filter = TRUE)
  soma_imputed <- impute(soma_pro, 0)
  DE <- do_limma_by_condition(soma_imputed, condition = "Sex", control_group = "F", treatment_group = "M")
  p <- plot_volcano(DE, label_col = "Genes")
  expect_s4_class(soma_pro, "SummarizedExperiment")
})
