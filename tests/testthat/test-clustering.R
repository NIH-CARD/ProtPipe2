dat <- data.table::fread("EXAMPLES/basic_example_data/iPSC.csv")
dat_pro <- create_se(dat)

test_that("correctly perform PCA clustering", {
  expect_error(plot_PCs(dat_pro))
  p <- plot_PCs(dat_pro %>% ProtPipe::impute_min())
  expect_s3_class(p, "ggplot")
})

test_that("correctly perform hierarchial clustering", {
  expect_error(plot_hierarchical_cluster(dat_pro))
  p <- plot_hierarchical_cluster(dat_pro %>% ProtPipe::impute_min())
  expect_s3_class(p, "ggplot")
})

test_that("correctly perform umap clustering", {
  expect_error(plot_umap(dat_pro))
  p <- plot_umap(dat_pro %>% ProtPipe::impute_min())
  expect_s3_class(p, "ggplot")
})
