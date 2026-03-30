test_that("PCA helpers return structured output on imputed data", {
  se <- load_basic_imputed_se()

  pcs <- ProtPipe::get_PCs(se, condition = "base_condition")

  expect_named(pcs, c("summary", "components"))
  expect_true(all(c("PC1", "PC2", "Sample", "Condition") %in% names(pcs$components)))

  p <- ProtPipe::plot_PCs(se, condition = "base_condition", pc_x = "PC1", pc_y = "PC2")
  expect_s3_class(p, "ggplot")
})

test_that("hierarchical clustering requires imputed data", {
  expect_error(
    ProtPipe::plot_hierarchical_cluster(
      load_basic_se(),
      dist_method = "euclidean",
      hclust_method = "complete"
    )
  )

  p <- ProtPipe::plot_hierarchical_cluster(
    load_basic_imputed_se(),
    dist_method = "euclidean",
    hclust_method = "complete"
  )
  expect_s3_class(p, "ggplot")
})

test_that("UMAP helpers return structured output on imputed data", {
  skip_if_not_installed("umap")

  set.seed(1)
  se <- load_basic_imputed_se()
  umap_df <- ProtPipe::get_umap(se, condition = "base_condition", neighbors = 5)

  expect_true(all(c("Sample", "UMAP1", "UMAP2", "Condition") %in% names(umap_df)))

  set.seed(1)
  p <- ProtPipe::plot_umap(se, condition = "base_condition", neighbors = 5)
  expect_s3_class(p, "ggplot")
})
