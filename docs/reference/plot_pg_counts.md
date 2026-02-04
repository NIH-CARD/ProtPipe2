# Plot Protein Group Counts

Generates a bar plot showing the number of identified protein groups per
sample. The plot can either show counts for each individual sample, or
show the summarized mean and standard deviation for samples grouped by a
condition.

## Usage

``` r
plot_pg_counts(object, condition = NULL)
```

## Arguments

- object:

  A `SummarizedExperiment` object.

- condition:

  Optional. A character string specifying a column name in the `colData`
  slot. If provided, samples will be grouped by this variable and the
  plot will show the mean and standard deviation of counts per group. If
  `NULL` (the default), each sample is plotted individually.

## Value

A `ggplot` object.
