# Create a SummarizedExperiment Object from Olink NPX Data

A wrapper function that processes a data frame of Olink NPX data,
optionally merges it with a sample condition/metadata file, and formats
the result into a SummarizedExperiment object for downstream analysis.

## Usage

``` r
create_se_from_olink(npx, condition = NULL, filter = TRUE)
```

## Arguments

- npx:

  A data frame containing Olink NPX data, typically from
  OlinkAnalyze::read_npx().

- condition:

  An optional data frame containing sample metadata. It is crucial that
  this file contains a column named 'SampleID' for merging.

- filter:

  A logical value (TRUE/FALSE) indicating whether to filter out samples
  with QC warnings. Defaults to FALSE.

## Value

A SummarizedExperiment object formatted for use with other package
functions.
