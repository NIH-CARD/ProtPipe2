# Calculate Principal Components Analysis (PCA)

Performs a principal component analysis (PCA) on the samples within a
`SummarizedExperiment` object.

## Usage

``` r
get_PCs(object, condition = NA)

# S4 method for class 'SummarizedExperiment'
get_PCs(object, condition = NA)
```

## Arguments

- object:

  A `SummarizedExperiment` object. The first assay is used.

- condition:

  A character string specifying the column name in `colData(object)` to
  use for grouping the samples in the output. If `NA` (the default),
  groups are inferred from sample names by removing numeric suffixes.

## Value

A list with two elements:

- `summary`:

  A data.table with the standard deviation, proportion of variance, and
  cumulative proportion for each component.

- `components`:

  A data frame with the first 5 PCA scores (PC1-PC5) for each sample,
  along with their assigned condition.

## Details

This function expects data that has already been cleaned and imputed.
Missing values (`NA`) will cause the PCA to fail. It is also highly
recommended to perform log-transformation and normalization before PCA.
The function will automatically remove proteins (rows) with zero
variance across samples before running the analysis.

## Functions

- `get_PCs(SummarizedExperiment)`: Method for SummarizedExperiment
  objects.

## Examples

``` r
# --- Create a sample SummarizedExperiment object ---
set.seed(123)
counts <- matrix(rnorm(100 * 6, mean = 10, sd = 2), nrow = 100, ncol = 6)
colnames(counts) <- paste0("Sample_", rep(c("A","B"), each=3), "_", 1:3)
sample_info <- data.frame(
  row.names = colnames(counts),
  Group = rep(c("A", "B"), each = 3)
)
se <- SummarizedExperiment(assays = list(log_intensities = counts),
                           colData = sample_info)
#> Error in SummarizedExperiment(assays = list(log_intensities = counts),     colData = sample_info): could not find function "SummarizedExperiment"

# --- Run PCA, specifying the group variable ---
pca_results <- get_PCs(se, condition = "Group")
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'get_PCs': object 'se' not found
head(pca_results$components)
#> Error: object 'pca_results' not found
```
