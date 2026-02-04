# Add a Processing Step to an Object's History

Appends a formatted log entry to the `processing` slot of an S4 object.
This function is intended to be used internally by preprocessing
functions (e.g., `log2_transform`, `impute_data`) to record their
actions.

## Usage

``` r
add_processing_step(object, log_entry)
```

## Arguments

- object:

  A `SummarizedExperiment` object with a processing_log entry in the
  metadata slot.

- log_entry:

  A named list containing metadata about the processing step. This list
  should ideally contain elements like `name`, `description`, and
  `parameters`.

## Value

The modified object with the `log_entry` appended to its
`processing_log` history.
