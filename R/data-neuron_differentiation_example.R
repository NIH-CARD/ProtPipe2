#' Bundled intensity table for the neuron differentiation example
#'
#' A data frame of protein group annotations and sample intensities from the
#' bundled neuron differentiation example dataset.
#'
#' @format A data frame with protein metadata columns followed by one column per
#'   sample intensity measurement.
#' @source `EXAMPLES/basic_example_data/iPSC.csv`
"neuron_differentiation_intensities"

#' Bundled sample metadata for the neuron differentiation example
#'
#' A data frame containing the sample identifiers and ordered differentiation
#' day labels used with `neuron_differentiation_intensities`.
#'
#' @format A data frame with 42 rows and 2 columns:
#' \describe{
#'   \item{SampleID}{Sample names matching the intensity columns in
#'   `neuron_differentiation_intensities`.}
#'   \item{differentiation_day}{An ordered factor with levels `day_0`, `day_3`,
#'   `day_7`, `day_10`, `day_14`, `day_21`, and `day_28`.}
#' }
#' @source `EXAMPLES/basic_example_data/iPSC.csv`
"neuron_differentiation_metadata"

#' Bundled iPSC stem cell marker genes
#'
#' A data frame of stem cell marker genes used in the README heatmap example.
#'
#' @format A data frame with one column, `Gene`.
#' @source `EXAMPLES/basic_example_data/stem_cell_gene.csv`
"ipsc_stem_cell_genes"
