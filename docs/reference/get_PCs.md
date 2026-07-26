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
se <- SummarizedExperiment::SummarizedExperiment(assays = list(log_intensities = counts),
                           colData = sample_info)

# --- Run PCA, specifying the group variable ---
pca_results <- get_PCs(se, condition = "Group")
head(pca_results$components)
#>                   PC1        PC2        PC3        PC4        PC5     Sample
#> Sample_A_1 -0.5406994 -6.9038248  1.0758652 -5.7618855  1.0603314 Sample_A_1
#> Sample_A_2  0.2890845 -3.8033696  0.7192578  7.2782128  2.9474868 Sample_A_2
#> Sample_A_3  2.3225423  1.7361782 -8.7966929 -0.7086532  0.1850881 Sample_A_3
#> Sample_B_1 -9.4317852  2.7544166  0.3131402  0.1892512 -1.6532474 Sample_B_1
#> Sample_B_2  3.0444035  5.9119874  4.0189782 -1.9767836  4.1974578 Sample_B_2
#> Sample_B_3  4.3164542  0.3046121  2.6694515  0.9798583 -6.7371168 Sample_B_3
#>            Condition
#> Sample_A_1         A
#> Sample_A_2         A
#> Sample_A_3         A
#> Sample_B_1         B
#> Sample_B_2         B
#> Sample_B_3         B
```
