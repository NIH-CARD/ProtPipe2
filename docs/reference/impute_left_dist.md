# Impute from a Down-Shifted Normal Distribution

Performs row-wise imputation by drawing random values from a normal
distribution that is shifted to the left and narrower than the
distribution of observed values.

This method assumes that missing values are primarily from proteins with
low abundance (i.e., below the detection limit). The default `shift` and
`scale` values are based on those used in the Perseus analysis platform.

## Usage

``` r
impute_left_dist(object, shift = 1.8, scale = 0.3)

# S4 method for class 'SummarizedExperiment'
impute_left_dist(object, shift = 1.8, scale = 0.3)
```

## Arguments

- object:

  A `SummarizedExperiment` object containing data with missing values.

- shift:

  A numeric value specifying how many standard deviations to shift the
  mean of the distribution for imputed values. Default is 1.8.

- scale:

  A numeric value to scale the standard deviation of the distribution
  for imputed values. Default is 0.3.

## Value

A `SummarizedExperiment` object with missing values imputed from a
simulated low-abundance distribution.

## Details

This imputation method should be applied to log-transformed data, as the
underlying assumption of a normal distribution is more appropriate in
log space. An exception is made for rows with only one or zero observed
values, where missing values are imputed with 0.

## Examples

``` r
# Create data with NAs, typically representing log-transformed values
raw_data <- data.frame(
  Gene = c("GENEA", "GENEB"),
  SampleA = c(25.1, 28.5),
  SampleB = c(25.5, 28.9),
  SampleC = c(NA, NA)
)

se <- create_se(raw_data)
#> `intensity_cols` not provided. Detecting numeric columns as intensity data.
#> Warning: `sample_metadata` not provided. Generating a basic version from column names.

# For reproducibility of the random imputation
set.seed(123)

imputed_obj <- impute_left_dist(se)
cat("Data after imputation:\n")
#> Data after imputation:
print(assay(imputed_obj))
#> Error in assay(imputed_obj): could not find function "assay"
```
