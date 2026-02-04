# Create a SummarizedExperiment Object from SomaScan Data

Processes SomaScan data from an ADAT object, optionally performs
filtering (e.g., on buffer controls), and formats the result into a
SummarizedExperiment object for downstream analysis.

## Usage

``` r
create_se_from_soma(adat, condition = NULL, filter = TRUE)
```

## Arguments

- adat:

  A SomaScan ADAT object, typically loaded using the SomaDataIO package.

- condition:

  An optional data frame containing sample metadata. This is passed to
  the downstream create_protdata function for integration.

- filter:

  A logical value (TRUE/FALSE) indicating whether to perform filtering
  on the data, including for buffer controls. Defaults to TRUE.

## Value

A SummarizedExperiment object formatted for use with other package
functions.
