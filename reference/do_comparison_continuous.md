# Perform limma differential expression for a continuous outcome

This function takes a SummarizedExperiment object, a numeric condition
column, and optional covariates, and performs limma-based differential
expression.

## Usage

``` r
do_comparison_continuous(object, condition, covariates = NULL)
```

## Arguments

- object:

  A SummarizedExperiment object

- condition:

  Column name in colData(object) used as the continuous predictor

- covariates:

  Optional covariates (must be in colData(object), cannot include
  `condition`)

## Value

A data frame with metadata, intensities, logFC, p-values
