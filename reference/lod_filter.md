# Filter SummarizedExperiment based on Limit of Detection (LOD)

This function filters a SummarizedExperiment object by setting values to
NA if they fall below a specified Limit of Detection (LOD). It searches
for the LOD values first in the rowData, and then in the assay columns
(samples).

## Usage

``` r
lod_filter(se, lod_col = "Buffer")
```

## Arguments

- se:

  A SummarizedExperiment object (e.g., proteomics data).

- lod_col:

  A character string specifying the name of the LOD column. Defaults to
  "Buffer".

## Value

The modified SummarizedExperiment object with values \< LOD set to NA.
