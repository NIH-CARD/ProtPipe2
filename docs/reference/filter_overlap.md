# Retain proteins present in a specified group of a `SummarizedExperiment` object

Filters the object to retain only proteins that are present (non-NA) in
at least one sample within each unique group of the specified condition.

## Usage

``` r
filter_overlap(object, condition_name)

# S4 method for class 'SummarizedExperiment,character'
filter_overlap(object, condition_name)
```

## Arguments

- object:

  A `SummarizedExperiment` object

- condition_name:

  A character string specifying the column name in the `colData` slot to
  group samples by.

## Value

A new, modified `SummarizedExperiment` object containing only the
overlapping proteins.
