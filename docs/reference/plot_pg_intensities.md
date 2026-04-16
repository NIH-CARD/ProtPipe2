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

# Generate the boxplot of intensities
p <- plot_pg_intensities(se)
if (interactive()) {
  print(p)
}
```
