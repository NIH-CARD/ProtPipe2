# Plot Boxplots of Sample Intensity Distributions

Generates boxplots to visualize and compare the distribution of protein
intensities for each sample. This is a common quality control plot to
check for normalization issues or identify potential outlier samples.

## Usage

``` r
plot_pg_intensities(object)

# S4 method for class 'SummarizedExperiment'
plot_pg_intensities(object)
```

## Arguments

- object:

  A `SummarizedExperiment` object.

## Value

A `ggplot` object showing a boxplot for each sample's log10-transformed
intensity distribution.

## Examples

``` r
# Create sample data
raw_data <- data.frame(
  Gene = c("GENEA", "GENEB", "GENEC", "GENED"),
  SampleA = c(100, 200, 150, 120),
  SampleB = c(110, 210, 160, 130), # Slightly higher than A
  SampleC = c(250, 500, 400, 300)   # Higher median and spread
)
se <- create_se(raw_data)
#> `intensity_cols` not provided. Detecting numeric columns as intensity data.
#> Warning: `sample_metadata` not provided. Generating a basic version from column names.
#> Error in SummarizedExperiment(assays = list(intensities = assay_data),     rowData = row_data, colData = col_data, metadata = list(creation_method = creation_method,         processing_log = list())): could not find function "SummarizedExperiment"

# Generate the boxplot of intensities
p <- plot_pg_intensities(se)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'plot_pg_intensities': object 'se' not found
if (interactive()) {
  print(p)
}
```
