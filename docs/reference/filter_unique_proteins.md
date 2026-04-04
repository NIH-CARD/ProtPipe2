# Removes duplicate analytes

Given a prot_meta column name, it will remove rows so that each value in
the column is unique. For each row the analyte with 1) the lowest
missing values and 2) the greatest median intensity will be kept.

## Usage

``` r
filter_unique_proteins(object, col = NULL)

# S4 method for class 'SummarizedExperiment,character'
filter_unique_proteins(object, col = NULL)
```

## Arguments

- object:

  A `SummarizedExperiment` object.

- col:

  a column name of the rowData slot

## Value

A `SummarizedExperiment` object.
