test_that("compare_protein returns a ggplot for SummarizedExperiment input", {
  p <- ProtPipe2::compare_protein(load_basic_se(), "U3KQP1")
  expect_s3_class(p, "ggplot")
})
