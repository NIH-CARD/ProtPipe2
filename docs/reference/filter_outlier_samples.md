# Filter Outlier Samples Based on Protein Counts

This function identifies and removes entire samples that are considered
outliers based on the total number of non-missing proteins identified
within them.

The method calculates the mean and standard deviation of protein counts
across all samples. A sample is flagged as an outlier if its total
protein count falls outside the range defined by a specified number of
standard deviations (`sds`) from the mean.

## Usage

``` r
filter_outlier_samples(object, sds = 3)

# S4 method for class 'SummarizedExperiment,numeric'
filter_outlier_samples(object, sds = 3)
```

## Arguments

- object:

  A `SummarizedExperiment` object.

- sds:

  The number of standard deviations from the mean protein count to use
  as the outlier threshold. Defaults to 3.

## Value

A `SummarizedExperiment` object with the outlier samples removed.
