# Calculate Pairwise Sample Correlations

Computes a pairwise correlation matrix for all samples based on their
protein intensity profiles. The result is returned in a long "tidy"
format data table.

## Usage

``` r
get_sample_correlation(object, method = "spearman")
```

## Arguments

- object:

  A `SummarizedExperiment` object.

- method:

  The correlation method to use. Can be one of `"spearman"` (the
  default), `"pearson"`, or `"kendall"`.

## Value

A `data.table` with three columns: `SampleA`, `SampleB`, and a third
column named for the method used (e.g., `Spearman`), containing the
pairwise correlation coefficients.
