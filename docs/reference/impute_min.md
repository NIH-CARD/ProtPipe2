# Impute Missing Values with the Row Minimum

Performs row-wise minimum imputation. For each protein (row), it
replaces missing values (`NA`, `NaN`) with the minimum observed value
found in that same row.

The imputation value can be scaled by a multiplicative factor `alpha`.

## Usage

``` r
impute_min(object, alpha = 1)

# S4 method for class 'SummarizedExperiment'
impute_min(object, alpha = 1)
```

## Arguments

- object:

  A `SummarizedExperiment` object containing data with missing values.

- alpha:

  A numeric scaling factor to multiply the row minimum by before
  imputation. Defaults to 1 (no scaling).

## Value

A `SummarizedExperiment` object with missing values imputed on a
per-protein basis.

## Examples

``` r
# Create data with different minimums and NAs in each row
raw_data <- data.frame(
  Gene = c("GENEA", "GENEB"),
  SampleA = c(100, 500),
  SampleB = c(200, 600),
  SampleC = c(NA, NA)
)

se <- create_se(raw_data)
#> `intensity_cols` not provided. Detecting numeric columns as intensity data.
#> Warning: `sample_metadata` not provided. Generating a basic version from column names.
cat("Original Data:\n")
#> Original Data:
print(assay(se))
#> Error in assay(se): could not find function "assay"

# Impute using the row minimum (alpha = 1)
# Row 1's NA becomes 100; Row 2's NA becomes 500.
imputed_obj <- impute_min(se)
cat("\nImputed with alpha = 1:\n")
#> 
#> Imputed with alpha = 1:
print(assay(imputed_obj))
#> Error in assay(imputed_obj): could not find function "assay"

# Impute using 90% of the row minimum
imputed_scaled <- impute_min(se, alpha = 0.9)
cat("\nImputed with alpha = 0.9:\n")
#> 
#> Imputed with alpha = 0.9:
print(assay(imputed_scaled))
#> Error in assay(imputed_scaled): could not find function "assay"
```
