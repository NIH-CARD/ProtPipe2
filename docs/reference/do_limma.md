# Perform limma differential expression on a ProtData object

This function takes a `SummarizedExperiment` object and two vectors
containing the column names of the treatment and control samples. It
filters out proteins found in fewer than 50% of samples, performs a
limma-based differential expression analysis, and returns a data frame
with metadata, intensities, log fold changes, and p-values.

## Usage

``` r
do_limma(object, treatment_samples, control_samples)
```

## Arguments

- object:

  A `SummarizedExperiment` object containing protein intensities,
  metadata, and condition info.

- treatment_samples:

  Character vector of column names representing treatment samples.

- control_samples:

  Character vector of column names representing control samples.

## Value

A data frame containing filtered proteins with metadata, intensity
values, log fold change, p-value, and adjusted p-value.
