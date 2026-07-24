# Plot a Hierarchical Clustering Dendrogram of Samples

Performs hierarchical clustering on the samples based on their protein
abundance profiles and generates a dendrogram plot using the `ggdendro`
package.

## Usage

``` r
plot_hierarchical_cluster(
  object,
  dist_method = "euclidean",
  hclust_method = "complete"
)

# S4 method for class 'SummarizedExperiment'
plot_hierarchical_cluster(
  object,
  dist_method = "euclidean",
  hclust_method = "complete"
)
```

## Arguments

- object:

  A `SummarizedExperiment` object. The data should be imputed.

- dist_method:

  The distance measure to be used by
  [`stats::dist`](https://rdrr.io/r/stats/dist.html). Common options
  include `"euclidean"`, `"maximum"`, `"manhattan"`. Defaults to
  `"euclidean"`.

- hclust_method:

  The agglomeration method to be used by
  [`stats::hclust`](https://rdrr.io/r/stats/hclust.html). Common options
  include `"complete"`, `"ward.D2"`, `"average"`. Defaults to
  `"complete"`.

## Value

A `ggplot` object representing the dendrogram, which can be further
customized.

## Details

This function expects clean, imputed data. Missing values (`NA`) will
cause an error. For meaningful biological results, it is highly
recommended to use data that has been log-transformed and normalized
before clustering.

## Functions

- `plot_hierarchical_cluster(SummarizedExperiment)`: Method for
  SummarizedExperiment objects.

## Examples

``` r
# Create a sample SummarizedExperiment object
raw_data <- data.frame(
  Gene = c("GENEA", "GENEB", "GENEC", "GENED"),
  SampleA = c(10, 20, 15, 12),
  SampleB = c(11, 21, 16, 13), # Similar to A
  SampleC = c(25, 10, 30, 5),
  SampleD = c(26, 11, 31, 6)  # Similar to C
)
se <- create_se(raw_data)
#> `intensity_cols` not provided. Detecting numeric columns as intensity data.
#> Warning: `sample_metadata` not provided. Generating a basic version from column names.

# Run with default methods. We expect A/B and C/D to cluster together.
p1 <- plot_hierarchical_cluster(se)
#> Error in assay(object): could not find function "assay"
if (interactive()) {
  print(p1)
}

# Run with different methods
p2 <- plot_hierarchical_cluster(se, dist_method = "manhattan", hclust_method = "ward.D2")
#> Error in assay(object): could not find function "assay"
if (interactive()) {
  print(p2)
}
```
