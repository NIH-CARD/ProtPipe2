# Filter Assay Based on Limit of Detection (LOD)

Filters the assay data of a SummarizedExperiment by comparing
intensities against a specific Limit of Detection (LOD). The LOD values
can be sourced either from a metadata column (rowData) or a specific
reference sample (e.g., a Buffer column in the assay).

## Usage

``` r
lod_filter(se, lod_col = "Buffer")

# S4 method for class 'SummarizedExperiment'
lod_filter(se, lod_col = "Buffer")
```

## Arguments

- se:

  A `SummarizedExperiment` object.

- lod_col:

  A character string indicating where to find the LOD values. The
  function first searches `rowData(se)`, then searches the column names
  of the assay (e.g., a specific buffer sample). Defaults to "Buffer".

## Value

A `SummarizedExperiment` object where values below the defined LOD are
replaced with `NA`.

## Functions

- `lod_filter(SummarizedExperiment)`: Method for SummarizedExperiment
  objects
