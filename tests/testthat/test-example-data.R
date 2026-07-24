test_that("protpipe_example_se loads as a SummarizedExperiment", {
  data("protpipe_example_se", package = "ProtPipe")

  expect_s4_class(protpipe_example_se, "SummarizedExperiment")
  expect_true("differentiation_day" %in% colnames(as.data.frame(SummarizedExperiment::colData(protpipe_example_se))))
  expect_equal(colnames(as.data.frame(SummarizedExperiment::colData(protpipe_example_se))), "differentiation_day")
})

test_that("bundled neuron differentiation tables load with expected structure", {
  data("neuron_differentiation_intensities", package = "ProtPipe")
  data("neuron_differentiation_metadata", package = "ProtPipe")
  data("ipsc_stem_cell_genes", package = "ProtPipe")

  expect_s3_class(neuron_differentiation_intensities, "data.frame")
  expect_s3_class(neuron_differentiation_metadata, "data.frame")
  expect_s3_class(ipsc_stem_cell_genes, "data.frame")

  expect_true(all(c("PG.ProteinGroups", "PG.Genes") %in% colnames(neuron_differentiation_intensities)))
  expect_identical(colnames(neuron_differentiation_metadata), c("SampleID", "differentiation_day"))
  expect_true(is.ordered(neuron_differentiation_metadata$differentiation_day))
  expect_identical(
    levels(neuron_differentiation_metadata$differentiation_day),
    c("day_0", "day_3", "day_7", "day_10", "day_14", "day_21", "day_28")
  )
  expect_identical(colnames(ipsc_stem_cell_genes), "Gene")
})

test_that("bundled neuron differentiation tables reconstruct a SummarizedExperiment", {
  data("neuron_differentiation_intensities", package = "ProtPipe")
  data("neuron_differentiation_metadata", package = "ProtPipe")

  se <- ProtPipe2::create_se(
    data = neuron_differentiation_intensities,
    sample_metadata = neuron_differentiation_metadata,
    creation_method = "bundled neuron differentiation example"
  )

  expect_s4_class(se, "SummarizedExperiment")
  expect_true("differentiation_day" %in% colnames(as.data.frame(SummarizedExperiment::colData(se))))
  expect_equal(ncol(SummarizedExperiment::assay(se)), nrow(neuron_differentiation_metadata))
})
