# Perform limma differential expression for a continuous outcome

This function takes a SummarizedExperiment object, a numeric condition
column, and optional covariates, and performs limma-based differential
expression.

## Usage

``` r
do_comparison_continuous(object, condition)

# S4 method for class 'SummarizedExperiment'
do_comparison_continuous(object, condition)
```

## Arguments

- object:

  A SummarizedExperiment object

- condition:

  Column name in colData(object) used as the continuous predictor

## Value

A data frame with metadata, intensities, logFC, p-values

## Methods (by class)

- `do_comparison_continuous(SummarizedExperiment)`: Method for
  SummarizedExperiment objects.
