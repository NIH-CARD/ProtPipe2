testthat::skip_if_not_installed("OlinkAnalyze")

# OlinkAnalyze ships npx_data1 itself; it is identical to the NPX export the
# suite used to read from EXAMPLES/olink/npx_data1.csv. The manifest is
# ProtPipe-specific, so it lives in inst/extdata.
npx <- OlinkAnalyze::npx_data1
meta <- read.delim(extdata_path("manifest.tsv"), sep = "\t")
out <- olink_all_output(npx)

test_that("correctly make a prot_data object from Olink without LOD filtering", {
  dat_pro <- create_se_from_olink(npx, meta, filter = T)
  t <- ProtPipe2::get_sample_correlation(dat_pro)
  tt <- ProtPipe2::plot_correlation_heatmap(dat_pro)
  expect_s4_class(dat_pro, "SummarizedExperiment")
  cols = ncol(SummarizedExperiment::assay(dat_pro))
  expect_equal(cols, 158)
  rows = nrow(SummarizedExperiment::assay(dat_pro))
  expect_equal(rows, 184)
})

test_that("correctly make a prot_data object from Olink with LOD filtering", {
  dat_pro <- create_se_from_olink(npx, meta, filter = TRUE)
  expect_s4_class(dat_pro, "SummarizedExperiment")
  cols = ncol(SummarizedExperiment::assay(dat_pro))
  expect_equal(cols, 158)
  rows = nrow(SummarizedExperiment::assay(dat_pro))
  expect_equal(rows, 184)
})

sampleIDs <- unique(npx$SampleID)
num_samples <- length(sampleIDs)
set.seed(42)
batches <- sample(c("Batch1", "Batch2", "Batch3"), num_samples, replace = TRUE)
ages <- sample(25:70, num_samples, replace = TRUE)
conditions <- sample(c("Control", "TreatmentA", "TreatmentB"), num_samples, replace = TRUE)

example_df <- data.frame(
  SampleID = sampleIDs,
  Batch = batches,
  Age = ages,
  Condition = conditions
)

output_filename <- tempfile(fileext = ".tsv")
write.table(
  example_df,
  file = output_filename,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Print a confirmation message
cat(paste0("\nData frame successfully saved to '", output_filename, "'\n"))
