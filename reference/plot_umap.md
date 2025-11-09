# Plot UMAP Dimensionality Reduction Results

Generates a 2D scatter plot of UMAP results to visualize sample
relationships. This is a convenient wrapper function that first calls
`get_umap` to calculate the coordinates and then plots them using
`ggplot2`.

## Usage

``` r
plot_umap(object, condition = NA, neighbors = 15, ...)

# S4 method for class 'SummarizedExperiment'
plot_umap(object, condition = NA, neighbors = 15, ...)
```

## Arguments

- object:

  A `SummarizedExperiment` object. The data should be imputed.

- condition:

  A character string specifying the column name in the `condition` slot
  to use for coloring the points. If `NULL` (the default), `get_umap`
  will attempt to guess groups from sample names.

- neighbors:

  The size of the local neighborhood for the UMAP calculation.

- ...:

  Additional arguments passed on to `get_umap`, and in turn to the
  [`umap::umap`](https://rdrr.io/pkg/umap/man/umap.html) function (e.g.,
  `min_dist`, `metric`).

## Value

A `ggplot` object, which can be further customized or printed.

## Details

This function requires imputed, log-transformed, and normalized data for
the best results, as these are the best inputs for the wrapped
`get_umap` function.

**Important:** UMAP is a stochastic algorithm. For a reproducible plot,
you **must set a seed** (e.g., `set.seed(123)`) *before* calling this
function.

This function requires the `umap` and `ggplot2` packages.

## Functions

- `plot_umap(SummarizedExperiment)`: Method for SummarizedExperiment
  objects.

## Examples

``` r
# This example requires the 'umap' package
if (requireNamespace("umap", quietly = TRUE)) {
  # Create a sample ProtData object
  raw_data <- data.frame(
    Gene = paste0("GENE", 1:10),
    Control_1 = rnorm(10, 10), Control_2 = rnorm(10, 10),
    Treat_A_1 = rnorm(10, 12), Treat_A_2 = rnorm(10, 12),
    Treat_B_1 = rnorm(10, 15), Treat_B_2 = rnorm(10, 15)
  )
  pd_obj <- create_protdata(dat = raw_data)

  # For reproducible UMAP results, set a seed!
  set.seed(42)

  # Generate the UMAP plot
  p <- plot_umap(pd_obj, n_neighbors = 3)

  # The plot can be printed in an interactive session
  if (interactive()) {
    print(p)
  }
}
#> Error: unable to find an inherited method for function ‘plot_umap’ for signature ‘object = "ProtData"’
```
