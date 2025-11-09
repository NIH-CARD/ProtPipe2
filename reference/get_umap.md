# Calculate UMAP Dimensionality Reduction

Performs Uniform Manifold Approximation and Projection (UMAP) on the
samples to generate a 2D embedding. This is a powerful non-linear method
for visualizing sample relationships. This function serves as a wrapper
for the [`umap::umap`](https://rdrr.io/pkg/umap/man/umap.html) function.

## Usage

``` r
get_umap(object, condition = NA, neighbors = 15, ...)

# S4 method for class 'SummarizedExperiment'
get_umap(object, condition = NA, neighbors = 15, ...)
```

## Arguments

- object:

  A `SummarizedExperiment` object. The data should be imputed.

- condition:

  A character string specifying the column name in the `condition` slot
  to use for labeling points. If `NULL` (the default), it will attempt
  to guess groups from sample names.

- neighbors:

  The size of the local neighborhood UMAP will look at. This is a key
  hyperparameter affecting the balance between local and global
  structure. Defaults to 15.

- ...:

  Additional arguments passed on to the
  [`umap::umap`](https://rdrr.io/pkg/umap/man/umap.html) function (e.g.,
  `min_dist`, `n_epochs`, `metric`).

## Value

A `data.table` with columns for `Sample`, `UMAP1`, `UMAP2`, and
`Condition`.

## Details

This function expects clean, imputed data. Missing values (`NA`) will
cause an error. For meaningful results, it is highly recommended to use
data that has been log-transformed and normalized.

**Important:** UMAP is a stochastic algorithm, meaning it will produce
slightly different results each time it is run. For reproducible
results, you **must set a seed** (e.g., `set.seed(123)`) before calling
this function.

This function requires the `umap` package to be installed from CRAN.

## Functions

- `get_umap(SummarizedExperiment)`: Method for SummarizedExperiment
  objects.
