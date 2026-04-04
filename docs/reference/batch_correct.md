# Correct for Batch Effects

This function uses the `removeBatchEffect` method from the limma package
to adjust for known batch effects in the data.

## Usage

``` r
batch_correct(object, batch_variable, bio_variables = NULL)

# S4 method for class 'SummarizedExperiment,character'
batch_correct(object, batch_variable, bio_variables = NULL)
```

## Arguments

- object:

  A `SummarizedExperiment` object. The first assay will be used.

- batch_variable:

  A single character string specifying the column name in
  `colData(object)` that contains the batch information.

- bio_variables:

  An optional character vector of column names in `colData(object)` that
  correspond to biological variables of interest. These effects will be
  preserved during the correction. If `NULL` (the default), only an
  intercept is protected.

## Value

A `SummarizedExperiment` object with a new assay named `"corrected"`
containing the batch-corrected data.

## Examples

``` r
# --- Create a sample SummarizedExperiment object ---
set.seed(123)
counts <- matrix(rnorm(100 * 6, mean = 10, sd = 2), nrow = 100, ncol = 6)
colnames(counts) <- paste0("Sample", 1:6)

sample_info <- data.frame(
  row.names = colnames(counts),
  Condition = rep(c("A", "B"), each = 3),
  Batch = rep(c("Batch1", "Batch2"), times = 3)
)

se <- SummarizedExperiment(assays = list(log_intensities = counts),
                           colData = sample_info)
#> Error in SummarizedExperiment(assays = list(log_intensities = counts),     colData = sample_info): could not find function "SummarizedExperiment"

# --- Run batch correction, preserving the "Condition" variable ---
corrected_se <- batch_correct(se,
                              batch_variable = "Batch",
                              bio_variables = "Condition")
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'batch_correct': object 'se' not found

# --- Inspect the result ---
assayNames(corrected_se) # Now includes "corrected"
#> Error in assayNames(corrected_se): could not find function "assayNames"
```
