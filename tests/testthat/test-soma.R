testthat::skip_if_not_installed("SomaDataIO")

# SomaDataIO ships its own example soma_adat, which exercises the same code
# paths as the 24 MB .adat the suite used to read from EXAMPLES/soma/.
dat <- SomaDataIO::example_data
data1 <- soma_all_output(dat)


test_that("create se object from somascan data without filtering", {
  soma_pro <- create_se_from_soma(dat, filter = FALSE)
  soma_imputed <- impute(soma_pro, 0)
  DE <- do_limma_binary(
    soma_imputed,
    condition = "Sex",
    control_group = "F",
    treatment_group = "M",
    covariates = NULL
  )
  p <- plot_volcano(DE, label_col = "Genes")
  expect_s4_class(soma_pro, "SummarizedExperiment")
})

test_that("create se object from somascan data with filtering", {
  soma_pro <- create_se_from_soma(dat, filter = TRUE)
  soma_imputed <- impute(soma_pro, 0)
  DE <- do_limma_binary(
    soma_imputed,
    condition = "Sex",
    control_group = "F",
    treatment_group = "M",
    covariates = NULL
  )
  p <- plot_volcano(DE, label_col = "Genes")
  expect_s4_class(soma_pro, "SummarizedExperiment")
})
