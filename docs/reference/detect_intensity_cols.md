# Identify Numeric Columns by Index

A helper function that scans a data frame and returns the integer
indices of the columns that contain numeric data.

This is primarily used internally by data loading functions to
automatically distinguish quantitative intensity columns from metadata
columns.

## Usage

``` r
detect_intensity_cols(df)
```

## Arguments

- df:

  A `data.frame` to be scanned.

## Value

An integer vector containing the column indices of all numeric columns
found in the input `data.frame`.
