# Mean Normalization of Proteomics Data

Corrects for sample-specific biases by scaling each sample (column) to
have the same mean abundance.

This method adjusts each sample's intensities by a multiplicative
factor, aligning the mean of every sample to a global mean calculated
from the entire dataset.

## Usage

``` r
mean_normalize(object)

# S4 method for class 'SummarizedExperiment'
mean_normalize(object)
```

## Arguments

- object:

  A `SummarizedExperiment` object

## Value

A `SummarizedExperiment` object with its data normalized by the
mean-centering method.

## Details

This multiplicative adjustment is typically used for raw,
non-log-transformed intensity data.
