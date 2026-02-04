# Get Protein Group Counts Per Sample

Calculates the number of non-missing protein groups (or features) for
each sample in a ProtData object.

## Usage

``` r
get_pg_counts(object)
```

## Arguments

- object:

  A `SummarizedExperiment` object.

## Value

A data frame with two columns: `Sample` (the sample names) and
`Protein_Groups` (the count of non-NA proteins for that sample).
