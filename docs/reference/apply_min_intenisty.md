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
se <- SummarizedExperiment(assays = list(counts = counts))
#> Error in SummarizedExperiment(assays = list(counts = counts)): could not find function "SummarizedExperiment"

# Apply LOD of 6
se_filtered <- applyLOD(se, lod = 6)
#> Error in applyLOD(se, lod = 6): could not find function "applyLOD"
assay(se_filtered)
#> Error in assay(se_filtered): could not find function "assay"
```
