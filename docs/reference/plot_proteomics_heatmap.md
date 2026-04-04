# Plot a Proteomics Heatmap

Generates a heatmap of protein expression. The function can operate in
two modes depending on the `condition` argument.

1.  **Individual Samples (default):** If `condition` is `NULL`, a
    heatmap of all individual samples is generated, with samples
    clustered by similarity.

2.  **Summarized Conditions:** If a `condition` is provided, the
    function first calculates the mean expression for each protein
    across the replicates within each condition group, then generates a
    heatmap of these mean values.

In both modes, the data is row-wise Z-scored before plotting.

## Usage

``` r
plot_proteomics_heatmap(
  object,
  protmeta_col = NULL,
  genes = NULL,
  title = NULL,
  condition = NULL,
  cluster_cols = FALSE,
  cluster_rows = FALSE
)

# S4 method for class 'SummarizedExperiment'
plot_proteomics_heatmap(
  object,
  protmeta_col = NULL,
  genes = NULL,
  title = NULL,
  condition = NULL,
  cluster_cols = FALSE,
  cluster_rows = FALSE
)
```

## Arguments

- object:

  A `SummarizedExperiment` object. It is recommended to use data that
  has been log-transformed and imputed.

- protmeta_col:

  A character string specifying the name of the column in the `rowData`
  slot to use for protein labels on the heatmap rows. Defaults to the
  first column of `rowData`.

- genes:

  Optional. A character vector of gene or protein names to subset the
  data to. Only these proteins will be included in the heatmap. The
  match is case-insensitive.

- title:

  Optional. A character string for the plot title.

- condition:

  Optional. A character string specifying a column name in the `colData`
  slot. If provided, triggers the summarized heatmap mode.

- cluster_cols:

  Logical indicating whether heatmap columns should be clustered.

- cluster_rows:

  Logical indicating whether heatmap rows should be clustered.

## Value

A `ggplot` object representing the heatmap.

## Functions

- `plot_proteomics_heatmap(SummarizedExperiment)`: Method for
  SummarizedExperiment objects.
