# Impute Missing Values with the Row Minimum

Performs row-wise minimum imputation. For each protein (row), it
replaces missing values (`NA`, `NaN`) with the minimum observed value
found in that same row.

The imputation value can be scaled by a multiplicative factor `alpha`.

## Usage

``` r
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

pd_obj <- create_protdata(dat = raw_data)
cat("Original Data:\n")
#> Original Data:
print(pd_obj@data)
#>   SampleA SampleB
#> 1     100     200
#> 2     500     600

# Impute using the row minimum (alpha = 1)
# Row 1's NA becomes 100; Row 2's NA becomes 500.
imputed_obj <- impute_min(pd_obj)
#> Error: unable to find an inherited method for function ‘impute_min’ for signature ‘object = "ProtData"’
cat("\nImputed with alpha = 1:\n")
#> 
#> Imputed with alpha = 1:
print(imputed_obj@data)
#> Error: object 'imputed_obj' not found

# Impute using 90% of the row minimum
imputed_scaled <- impute_min(pd_obj, alpha = 0.9)
#> Error: unable to find an inherited method for function ‘impute_min’ for signature ‘object = "ProtData"’
cat("\nImputed with alpha = 0.9:\n")
#> 
#> Imputed with alpha = 0.9:
print(imputed_scaled@data)
#> Error: object 'imputed_scaled' not found
```
