# Calculate Coefficient of Variation (CV) for Protein Groups

Calculates the coefficient of variation (CV) for each protein, grouped
by an experimental condition. The CV is a standardized measure of
dispersion, calculated as the standard deviation divided by the mean,
and is often used to assess the reproducibility of replicates within a
condition group.

## Usage

``` r
get_CVs(object, condition, min_samples = 2)

# S4 method for class 'SummarizedExperiment'
get_CVs(object, condition, min_samples = 2)
```

## Arguments

- object:

  A `SummarizedExperiment` object.

- condition:

  A character string specifying the column name in the `colData` slot
  which defines the replicate groups.

- min_samples:

  The minimum number of replicates in a group required to calculate CVs.
  Groups with fewer samples will be skipped. Defaults to 2.

## Value

A data frame in long format with columns for `Protein`, `CV` (in
percent), and `Condition`.
