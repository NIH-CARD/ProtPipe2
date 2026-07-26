# Impute Missing Values with the Row Minimum

Performs row-wise minimum imputation. For each protein (row), it
replaces missing values (`NA`, `NaN`) with the minimum observed value
found in that same row.

The imputation value can be scaled by a multiplicative factor `alpha`.

Negative values are treated as missing and converted to `NA` before
imputation. Proteins with no observed values in any sample (all-NA rows)
are imputed with `0`.

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
print(SummarizedExperiment::assay(se))
#>      SampleA SampleB
#> [1,]     100     200
#> [2,]     500     600

# Impute using the row minimum (alpha = 1)
# Row 1's NA becomes 100; Row 2's NA becomes 500.
imputed_obj <- impute_min(se)
cat("\nImputed with alpha = 1:\n")
#> 
#> Imputed with alpha = 1:
print(SummarizedExperiment::assay(imputed_obj))
#>      SampleA SampleB
#> [1,]     100     200
#> [2,]     500     600

# Impute using 90% of the row minimum
imputed_scaled <- impute_min(se, alpha = 0.9)
cat("\nImputed with alpha = 0.9:\n")
#> 
#> Imputed with alpha = 0.9:
print(SummarizedExperiment::assay(imputed_scaled))
#>      SampleA SampleB
#> [1,]     100     200
#> [2,]     500     600
```
