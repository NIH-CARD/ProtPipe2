# Safely Convert Character Columns to Numeric Type

This function inspects each column of a data frame. If a character
column contains only values that can be interpreted as numbers
(including integers, decimals, scientific notation, and the strings "NA"
or "NaN"), it converts that entire column to the numeric type. Columns
containing non-numeric text are left unchanged.

## Usage

``` r
convert_numeric_cols(df)
```

## Arguments

- df:

  A data frame whose character columns will be checked for potential
  conversion to numeric.

## Value

The input data frame, with any eligible character columns converted to
the numeric class.
