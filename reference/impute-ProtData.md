# Impute Missing Values with a Constant

Replaces all missing values (`NA` and `NaN`) in the numeric columns of
the data with a single, user-specified constant.

## Usage

``` r
impute(object, value)
```

## Arguments

- object:

  A `SummarizedExperiment` object containing data with missing values.

- value:

  The numeric constant to use for replacing `NA` and `NaN` values.

## Value

A `SummarizedExperiment` object with missing values imputed.

## Examples

``` r
# Create a sample data frame with metadata and a missing numeric value
raw_data <- data.frame(
  Protein.ID = c("P02768", "P01023", "P60709"),
  Gene.Name = c("ALB", "A2M", "ACTB"),
  Sample_A = c(1.2e6, 2.3e6, NA),
  Sample_B = c(1.4e6, 2.6e6, 4.8e6)
)

# Use the constructor to create a ProtData object.
# The constructor will automatically separate metadata from numeric data.
pd_obj <- create_protdata(dat = raw_data)

# Impute the NA with 0
imputed_obj <- impute(pd_obj, value = 0)
#> Error: unable to find an inherited method for function ‘impute’ for signature ‘object = "ProtData", value = "numeric"’

# View the imputed data slot
print(imputed_obj@data)
#> Error: object 'imputed_obj' not found
```
