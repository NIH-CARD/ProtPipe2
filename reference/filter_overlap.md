# get_overlap method for protdata class

Filters the protdata object to retain only proteins that are present
(non-NA) in at least one sample within each unique group of the
specified condition.

## Usage

``` r
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
