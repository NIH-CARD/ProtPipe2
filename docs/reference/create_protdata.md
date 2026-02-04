# Create a ProtData Object

This function creates an instance of the ProtData class.

## Usage

``` r
create_protdata(
  dat,
  intensity_cols = NULL,
  condition = NULL,
  method = "Unknown"
)
```

## Arguments

- dat:

  A data frame containing proteomics data (proteins are rows, samples
  are columns).

- condition:

  A data frame containing conditions of the samples. Rownames should
  match colnames of data. Optional.

- method:

  A character string describing the method used for generating the data.
  Optional.

## Value

An instance of the ProtData class.
