# Perform t-test differential expression on a SummarizedExperiment

This function takes a `SummarizedExperiment` object and uses group
labels in `colData(object)` to perform differential expression between a
treatment group and a control group using per-protein t-tests. Optional
covariates can be regressed out before testing.

## Usage

``` r
do_t_test_binary(
  object,
  condition,
  treatment_group,
  control_group,
  covariates = NULL
)

# S4 method for class 'SummarizedExperiment'
do_t_test_binary(
  object,
  condition,
  treatment_group,
  control_group,
  covariates = NULL
)
```

## Arguments

- object:

  A `SummarizedExperiment` object containing protein intensities and
  metadata.

- condition:

  String: column name in `colData(object)` that holds group labels.

- treatment_group:

  String: name of treatment group in `condition`.

- control_group:

  String: name of control group in `condition`.

- covariates:

  Optional character vector: column names in `colData(object)` to
  regress out before testing.

## Value

A data frame with metadata, intensities, log fold change, p-values, and
adjusted p-values.

## Functions

- `do_t_test_binary(SummarizedExperiment)`: Method for
  SummarizedExperiment objects.
