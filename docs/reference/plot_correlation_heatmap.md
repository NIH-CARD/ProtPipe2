# Plot a Sample Correlation Heatmap

Creates a heatmap of pairwise sample correlations, calculated via the
Spearman method. The function includes a "smart sorting" feature to
correctly order axes based on numeric parts of labels (e.g., "Day1",
"Day2", "Day10").

## Usage

``` r
plot_correlation_heatmap(object, order_by = NULL, label_by = NULL)

# S4 method for class 'SummarizedExperiment'
plot_correlation_heatmap(object, order_by = NULL, label_by = NULL)
```

## Arguments

- object:

  A `SummarizedExperiment` object.

- order_by:

  Optional. A character string specifying a column in the condition slot
  by which to order the samples on the heatmap axes. If the column
  contains numeric parts, a natural sort is applied.

- label_by:

  Optional. A character string specifying a column in the condition slot
  to use for relabeling the heatmap axes.

## Value

A `ggplot` object representing the heatmap.
