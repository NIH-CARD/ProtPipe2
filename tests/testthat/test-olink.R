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

test_that("Olink NPX is not log2 transformed again by downstream analysis", {
  dat_pro <- create_se_from_olink(npx, meta, filter = FALSE)

  # NPX already is a log2 scale, so import must record the step...
  expect_true(ProtPipe2::has_step(dat_pro, "log2_transform"))

  # ...and the negative NPX values that log2(x + 1) would turn into NaN
  # must survive untouched.
  before <- SummarizedExperiment::assay(dat_pro)
  expect_true(any(before < 0, na.rm = TRUE))
  expect_false(any(is.nan(before)))

  imputed <- ProtPipe2::impute(dat_pro, 0)
  DE <- ProtPipe2::do_limma_binary(
    imputed,
    condition = "Condition",
    control_group = "Control",
    treatment_group = "TreatmentA",
    covariates = NULL
  )

  sample_col <- intersect(colnames(DE), colnames(before))[1]
  expect_equal(
    sort(DE[[sample_col]]),
    sort(as.numeric(SummarizedExperiment::assay(imputed)[, sample_col]))
  )
})

test_that("SomaScan RFU is still log2 transformed (linear scale)", {
  skip_if_not_installed("SomaDataIO")
  soma <- create_se_from_soma(SomaDataIO::example_data, filter = FALSE)
  expect_false(ProtPipe2::has_step(soma, "log2_transform"))
})

test_that("Olink objects are flagged log2 regardless of import path", {
  # The Shiny app does not call create_se_from_olink(); it pipes read_NPX()
  # through olink_all_output() into create_se() directly, so the flag has to
  # come from creation_method rather than from the wrapper.
  wide <- ProtPipe2::olink_all_output(npx)$data
  via_app <- ProtPipe2::create_se(
    dat = wide,
    intensity_cols = 5:ncol(wide),
    creation_method = "Olink"
  )
  expect_true(ProtPipe2::has_step(via_app, "log2_transform"))
  expect_length(S4Vectors::metadata(via_app)$processing_log, 1)

  via_wrapper <- create_se_from_olink(npx, meta, filter = FALSE)
  expect_true(ProtPipe2::has_step(via_wrapper, "log2_transform"))
  # the wrapper must not record the step a second time
  expect_length(S4Vectors::metadata(via_wrapper)$processing_log, 1)

  # non-Olink platforms are untouched
  for (m in c("SomaScan", "Standard Matrix", "Unknown")) {
    se <- ProtPipe2::create_se(dat = wide, intensity_cols = 5:ncol(wide), creation_method = m)
    expect_false(ProtPipe2::has_step(se, "log2_transform"))
  }
})
