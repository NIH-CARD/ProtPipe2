# Plot Coefficient of Variation (CV) Distributions

Generates a violin or jitter plot to visualize the distribution of
protein CVs for each experimental condition. This plot is useful for
comparing the reproducibility and technical variance across different
sample groups.

## Usage

``` r
plot_CVs(object, condition, plot_type = "violin")

# S4 method for class 'SummarizedExperiment'
plot_CVs(object, condition, plot_type = "violin")
```

## Arguments

- object:

  A `SummarizedExperiment` object.

- condition:

  A character string specifying the column name in the `colData` slot
  which defines the replicate groups.

- plot_type:

  The type of plot to generate. Must be either `"violin"` (the default,
  which also includes a narrow boxplot) or `"jitter"`.

## Value

A `ggplot` object.
