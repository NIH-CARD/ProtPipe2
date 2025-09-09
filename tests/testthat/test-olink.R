library(testthat)
setwd("../../..")

load("EXAMPLES/olink/npx_data1.rda")
file = "EXAMPLES/olink/npx_data1.csv"
metafile <- "EXAMPLES/olink/manifest.tsv"

npx <- OlinkAnalyze::read_NPX(file)

meta <- read.delim(metafile, sep = "\t")

test_that("correctly make a prot_data object from Olink without LOD filtering", {
  dat_pro <- create_protdata_from_olink(npx, meta, filter = T)
  t <- ProtPipe::get_spearman(dat_pro)
  tt <- ProtPipe::plot_correlation_heatmap(dat_pro)
  expect_s4_class(dat_pro, "ProtData")
  cols = ncol(dat_pro@data)
  expect_equal(cols, 158)
  rows = nrow(dat_pro@data)
  expect_equal(rows, 184)
})

test_that("correctly make a prot_data object from Olink with LOD filtering", {
  dat_pro <- create_protdata_from_olink(npx)
  expect_s4_class(dat_pro, "ProtData")
  cols = ncol(dat_pro@data)
  expect_equal(cols, 158)
  rows = nrow(dat_pro@data)
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

output_filename <- "EXAMPLES/olink/manifest.tsv"
write.table(
  example_df,
  file = output_filename,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Print a confirmation message
cat(paste0("\nData frame successfully saved to '", output_filename, "'\n"))
