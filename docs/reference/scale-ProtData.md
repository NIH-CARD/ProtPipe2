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
# Create sample data with proteins as rows
df <- data.frame(
  SampleA = c(100, 250, 50),
  SampleB = c(120, 200, 100),
  SampleC = c(110, 225, 75),
  row.names = c("Protein1", "Protein2", "Protein3")
)

conditions <- data.frame(
  row.names = colnames(df),
  group = c("Control", "Treatment", "Control")
)

prot_obj <- new("ProtData",
                data = df,
                condition = conditions,
                method = "MS")

# Check the means of each protein (row) before scaling
rowMeans(prot_obj@data)
#> Protein1 Protein2 Protein3 
#>      110      225       75 

# Apply the scaling method
scaled_prot_obj <- scale(prot_obj)
#> Error in as.vector(data): no method for coercing this S4 class to a vector

# The new data has row means near zero and row standard deviations of one
print(scaled_prot_obj@data)
#> Error: object 'scaled_prot_obj' not found
cat("Row means after scaling:\n")
#> Row means after scaling:
print(rowMeans(scaled_prot_obj@data))
#> Error: object 'scaled_prot_obj' not found
cat("\nRow standard deviations after scaling:\n")
#> 
#> Row standard deviations after scaling:
print(apply(scaled_prot_obj@data, 1, sd))
#> Error: object 'scaled_prot_obj' not found
```
