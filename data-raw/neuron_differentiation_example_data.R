# Build packaged example datasets from the bundled neuron differentiation CSV.
# Run from the project root.

raw_dat <- utils::read.csv(
  "EXAMPLES/basic_example_data/iPSC.csv",
  check.names = FALSE
)

sample_ids <- colnames(raw_dat)[vapply(raw_dat, is.numeric, logical(1))]
day_numbers <- as.integer(sub("^\\[[0-9]+\\] Day([0-9]+)_.*$", "\\1", sample_ids))

if (anyNA(day_numbers)) {
  stop("Failed to parse differentiation days from example sample names.")
}

day_levels <- paste0("day_", sort(unique(day_numbers)))

neuron_differentiation_intensities <- raw_dat
neuron_differentiation_metadata <- data.frame(
  SampleID = sample_ids,
  differentiation_day = factor(
    paste0("day_", day_numbers),
    levels = day_levels,
    ordered = TRUE
  ),
  stringsAsFactors = FALSE
)

ipsc_stem_cell_genes <- utils::read.csv(
  "EXAMPLES/basic_example_data/stem_cell_gene.csv",
  check.names = FALSE
)

if (!dir.exists("data")) {
  dir.create("data")
}

save(
  neuron_differentiation_intensities,
  file = "data/neuron_differentiation_intensities.rda",
  compress = "bzip2"
)
save(
  neuron_differentiation_metadata,
  file = "data/neuron_differentiation_metadata.rda",
  compress = "bzip2"
)
save(
  ipsc_stem_cell_genes,
  file = "data/ipsc_stem_cell_genes.rda",
  compress = "bzip2"
)
