# Plot Principal Component Analysis Results

Generates a scatter plot of principal components to visualize the
relationships between samples. This function serves as a convenient
wrapper that first calls `calculate_pca` and then plots the results
using `ggplot2`.

## Usage

``` r
plot_PCs(object, condition = NA, pc_x = "PC1", pc_y = "PC2")

# S4 method for class 'SummarizedExperiment'
plot_PCs(object, condition = NA, pc_x = "PC1", pc_y = "PC2")
```

## Arguments

- object:

  A `SummarizedExperiment` object.

- condition:

  A character string specifying the column name in the `condition` slot
  to use for coloring the points. If `NA` (the default), it will attempt
  to guess groups from sample names.

- pc_x:

  A character string for the principal component to plot on the x-axis
  (e.g., `"PC1"`). Defaults to `"PC1"`.

- pc_y:

  A character string for the principal component to plot on the y-axis
  (e.g., `"PC2"`). Defaults to `"PC2"`.

## Value

A `ggplot` object, which can be further customized or printed.

## Details

As this function calls `calculate_pca` internally, the data in the
object must be imputed first. It is also recommended to use
log-transformed and normalized data for the best results.

## Functions

- `plot_PCs(SummarizedExperiment)`: Method for SummarizedExperiment
  objects.

## Examples

``` r
# Create a sample ProtData object with missing data
raw_data <- data.frame(
  Gene = c("GENEA", "GENEB", "GENEC", "GENED", "GENEE", "GENEF"),
  Control_1 = c(10, 11, 12, 13, 14, 15),
  Control_2 = c(10.5, 11.5, 12.5, NA, 14.5, 15.5),
  Treatment_1 = c(15, 16, 17, 18, 19, 20),
  Treatment_2 = c(15.5, 16.5, 17.5, 18.5, 19.5, 20.5)
)
cond_df <- data.frame(
   SampleID = c("Control_1", "Control_2", "Treatment_1", "Treatment_2"),
   group = c("Control", "Control", "Treatment", "Treatment")
)
pd_obj <- create_protdata(dat = raw_data, condition = cond_df)

# Impute missing values before plotting
pd_obj_imputed <- impute(pd_obj, value = 13.5)
#> Error: unable to find an inherited method for function ‘impute’ for signature ‘object = "ProtData", value = "numeric"’

# Generate the plot of PC1 vs PC2
p1 <- plot_pca(pd_obj_imputed, condition = "group")
#> Error in plot_pca(pd_obj_imputed, condition = "group"): could not find function "plot_pca"
if (interactive()) {
  print(p1)
}

# Generate a plot of PC1 vs PC3
p2 <- plot_pca(pd_obj_imputed, condition = "group", pc_x = "PC1", pc_y = "PC3")
#> Error in plot_pca(pd_obj_imputed, condition = "group", pc_x = "PC1", pc_y = "PC3"): could not find function "plot_pca"
if (interactive()) {
  print(p2)
}
```
