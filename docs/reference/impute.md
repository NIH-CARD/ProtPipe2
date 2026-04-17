# Impute Missing Values with a Constant

Replaces all missing values (`NA` and `NaN`) in the numeric columns of
the data with a single, user-specified constant.

## Usage

``` r
impute(object, value)

# S4 method for class 'SummarizedExperiment,numeric'
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

# Create a SummarizedExperiment object
se <- create_se(raw_data)
#> `intensity_cols` not provided. Detecting numeric columns as intensity data.
#> Warning: `sample_metadata` not provided. Generating a basic version from column names.

# Impute the NA with 0
imputed_obj <- impute(se, value = 0)

# View the imputed assay data
print(assay(imputed_obj))
#> Error in assay(imputed_obj): could not find function "assay"
```
