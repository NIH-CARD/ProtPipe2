# Filter Proteins by Percentage of Valid Values

A generic function to filter proteins (rows) from a proteomics data
object based on the percentage of samples where they have valid
observations.

## Usage

``` r
filter_proteins_by_percent(object, percent)

# S4 method for class 'SummarizedExperiment,numeric'
filter_proteins_by_percent(object, percent)
```

## Arguments

- object:

  A `SummarizedExperiment` object.

- percent:

  The minimum percentage (from 0 to 100) of samples in which a protein
  must have a valid (non-NA) value to be retained.

## Value

A modified object of the same class with proteins filtered.
