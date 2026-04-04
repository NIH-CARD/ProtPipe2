# Z-Score Normalization for Proteins Across Samples

This method performs a Z-score transformation on a protein-wise
(row-wise) basis. For each protein, it calculates the mean and standard
deviation of its abundance across all samples and then scales the values
accordingly.

This is a standard transformation for visualizing expression patterns in
a heatmap, as it highlights the relative change of each protein across
samples, independent of its absolute abundance.

## Usage

``` r
z_score(object)

# S4 method for class 'SummarizedExperiment'
z_score(object)
```

## Arguments

- object:

  A `SummarizedExperiment` object

## Value

object A `SummarizedExperiment` object where the abundance values in the
`assay` slot have been replaced by their row-wise Z-scores.

## See also

[`base::scale()`](https://rdrr.io/r/base/scale.html)

## Examples

``` r
# Create sample data with proteins as rows and samples as columns
df <- data.frame(
  Protein = c("Protein1", "Protein2", "Protein3"),
  SampleA = c(100, 250, 50),
  SampleB = c(120, 200, 100),
  SampleC = c(110, 225, 75)
)

sample_info <- data.frame(
  SampleID = c("SampleA", "SampleB", "SampleC"),
  group = c("Control", "Treatment", "Control")
)

se <- create_se(df, sample_metadata = sample_info)
#> `intensity_cols` not provided. Detecting numeric columns as intensity data.
#> Error in SummarizedExperiment(assays = list(intensities = assay_data),     rowData = row_data, colData = col_data, metadata = list(creation_method = creation_method,         processing_log = list())): could not find function "SummarizedExperiment"

# Check the means of each protein (row) before scaling
rowMeans(assay(se))
#> Error in assay(se): could not find function "assay"

# Apply the scaling method
scaled_se <- z_score(se)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'z_score': object 'se' not found

# The new data has row means near zero and row standard deviations of one
print(assay(scaled_se))
#> Error in assay(scaled_se): could not find function "assay"
cat("Row means after scaling:\n")
#> Row means after scaling:
print(rowMeans(assay(scaled_se)))
#> Error in assay(scaled_se): could not find function "assay"
cat("\nRow standard deviations after scaling:\n")
#> 
#> Row standard deviations after scaling:
print(apply(assay(scaled_se), 1, sd))
#> Error in assay(scaled_se): could not find function "assay"
```
