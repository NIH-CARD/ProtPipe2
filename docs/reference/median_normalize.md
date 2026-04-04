# Median Normalization of Proteomics Data

Performs median normalization on the quantitative data within a
`SummarizedExperiment` object. This method corrects for systematic,
sample-specific biases (e.g., differences in sample loading or
instrument sensitivity) to make the samples more comparable.

The function operates by first calculating the median intensity for each
sample (column). It then calculates the median of these medians (the
"global median"). Finally, it scales the intensities in each sample by a
multiplicative factor so that every sample has the same median
abundance, equal to the global median.

## Usage

``` r
median_normalize(object)

# S4 method for class 'SummarizedExperiment'
median_normalize(object)
```

## Arguments

- object:

  A `SummarizedExperiment` object

## Value

A `SummarizedExperiment` object with the abundance data in the `assay`
slot normalized by the median-centering method.

## Details

This type of multiplicative adjustment is typically appropriate for raw,
non-log-transformed intensity data.
