# Apply Limit of Detection Threshold

Filters assay data by converting values below a specified Limit of
Detection (LOD) to `NA`.

## Usage

``` r
apply_min_intenisty(object, lod)

# S4 method for class 'SummarizedExperiment'
apply_min_intenisty(object, lod)
```

## Arguments

- object:

  A `SummarizedExperiment` object containing a single proteomics assay.

- lod:

  A single numeric value representing the Limit of Detection.
  Intensities below this value will be replaced with `NA`.

## Value

The original `object` with the assay matrix updated to contain `NA`
values where intensities were below the `lod`.

## Examples

``` r
# Create a mock SummarizedExperiment
counts <- matrix(c(10, 5, 2, 8, 1, 15), nrow = 3, ncol = 2)
se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts))

# Apply LOD of 6
se_filtered <- apply_min_intenisty(se, lod = 6)
SummarizedExperiment::assay(se_filtered)
#>      [,1] [,2]
#> [1,]   10    8
#> [2,]   NA   NA
#> [3,]   NA   15
```
