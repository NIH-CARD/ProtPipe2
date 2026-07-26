# Perform ANOVA-style differential expression across multiple groups

This function takes a `SummarizedExperiment` object and uses group
labels in `colData(object)` to test, for each protein, whether mean
abundance differs across *all* levels of `condition` at once. It uses
limma's moderated F-test (empirical Bayes-shrunk variance, as in
`do_limma_binary`) rather than a classical per-protein
[`aov()`](https://rdrr.io/r/stats/aov.html), which is more robust when
group sizes are small. `condition` must have at least 3 groups; for
exactly 2 groups, use `do_limma_binary` instead.

## Usage

``` r
do_anova(object, condition, covariates = NULL)

# S4 method for class 'SummarizedExperiment'
do_anova(object, condition, covariates = NULL)
```

## Arguments

- object:

  A `SummarizedExperiment` object containing protein intensities and
  metadata.

- condition:

  String: column name in `colData(object)` that holds group labels (\>=
  3 groups).

- covariates:

  Optional character vector: column names in `colData(object)` to use as
  covariates.

## Value

A data frame with metadata, intensities, one log fold change column per
non-reference group (versus an arbitrary baseline group), the omnibus
`F` statistic, `P.Value`, and `adj.P.Val`.

## Functions

- `do_anova(SummarizedExperiment)`: Method for SummarizedExperiment
  objects.
