test_that("protpipe_example_se loads as a SummarizedExperiment", {
  data("protpipe_example_se", package = "ProtPipe")

  expect_s4_class(protpipe_example_se, "SummarizedExperiment")
  expect_true("differentiation_day" %in% colnames(as.data.frame(SummarizedExperiment::colData(protpipe_example_se))))
  expect_equal(colnames(as.data.frame(SummarizedExperiment::colData(protpipe_example_se))), "differentiation_day")
})
