# Bundled sample metadata for the neuron differentiation example

A data frame containing the sample identifiers and ordered
differentiation day labels used with
`neuron_differentiation_intensities`.

## Usage

``` r
neuron_differentiation_metadata
```

## Format

A data frame with 42 rows and 2 columns:

- `SampleID`:

  Sample names matching the intensity columns in
  `neuron_differentiation_intensities`.

- `differentiation_day`:

  An ordered factor with levels `day_0`, `day_3`, `day_7`, `day_10`,
  `day_14`, `day_21`, and `day_28`.

## Source

`EXAMPLES/basic_example_data/iPSC.csv`
