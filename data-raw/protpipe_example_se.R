# Build the packaged SummarizedExperiment example from the bundled iPSC CSV.
# Run from the project root.

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload must be installed to build protpipe_example_se.")
}

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("data.table must be installed to build protpipe_example_se.")
}

library(SummarizedExperiment)

pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)

raw_dat <- data.table::fread("EXAMPLES/basic_example_data/iPSC.csv")

protpipe_example_se <- ProtPipe::create_se(
  data = raw_dat,
  creation_method = "basic_example_data/iPSC.csv"
)

colnames(SummarizedExperiment::colData(protpipe_example_se)) <- "differentiation_day"

if (!dir.exists("data")) {
  dir.create("data")
}

save(protpipe_example_se, file = "data/protpipe_example_se.rda", compress = "bzip2")
