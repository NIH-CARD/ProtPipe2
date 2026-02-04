# This function takes a data frame of proteomics data and its corresponding sample metadata to construct a SummarizedExperiment object. It handles detection of intensity columns, validation, and synchronization of metadata.

This function takes a data frame of proteomics data and its
corresponding sample metadata to construct a SummarizedExperiment
object. It handles detection of intensity columns, validation, and
synchronization of metadata.

## Usage

``` r
create_se(
  data,
  sample_metadata = NULL,
  intensity_cols = NULL,
  creation_method = "Unknown"
)
```

## Arguments

- data:

  A data frame containing both protein metadata and intensity values.

- sample_metadata:

  A data frame for sample metadata. It must contain a 'SampleID' column
  that matches the intensity column headers in 'data'. If NULL (the
  default), a basic metadata table is generated from the column names.

- intensity_cols:

  An optional character or numeric vector specifying which columns in
  'data' are the intensity/abundance columns. If NULL, the function will
  attempt to autodetect them as all numeric columns.

- creation_method:

  A character string to log how the data was created (e.g.,
  "MaxQuant_LFQ"). Stored in the object's metadata.

## Value

A `SummarizedExperiment` object.

## Examples

``` r
# --- Sample Data ---
test_data <- data.frame(
  ProteinID = c("P02768", "P01023", "P10636"),
  Gene = c("ALB", "A2M", "CALM1"),
  Sample_A_1 = c(10.1, 12.5, 13.2),
  Sample_A_2 = c(10.5, 12.8, 13.9),
  Sample_B_1 = c(15.2, 18.1, 19.0),
  Sample_B_2 = c(15.8, 18.5, 19.9)
)

sample_info <- data.frame(
  SampleID = c("Sample_A_1", "Sample_A_2", "Sample_B_1", "Sample_B_2"),
  Condition = c("A", "A", "B", "B")
)

# --- Function Call ---
se <- create_se(
  data = test_data,
  sample_metadata = sample_info
)
#> `intensity_cols` not provided. Detecting numeric columns as intensity data.
#> Error in SummarizedExperiment(assays = list(intensities = assay_data),     rowData = row_data, colData = col_data, metadata = list(creation_method = creation_method,         processing_log = list())): could not find function "SummarizedExperiment"

# --- Inspect the new object ---
print(se)
#> Error: object 'se' not found
assayNames(se)
#> Error in assayNames(se): could not find function "assayNames"
head(rowData(se))
#> Error in rowData(se): could not find function "rowData"
colData(se)
#> Error in colData(se): could not find function "colData"
metadata(se)
#> Error in metadata(se): could not find function "metadata"
```
