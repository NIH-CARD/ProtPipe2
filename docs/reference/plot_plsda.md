# Generate a PLS-DA scores plot for a ProtData object

Generate a PLS-DA scores plot for a ProtData object

## Usage

``` r
plot_plsda(object, group_variable, n_components = 2)

# S4 method for class 'ProtData,character'
plot_plsda(object, group_variable, n_components = 2)
```

## Arguments

- object:

  The ProtData object.

- group_variable:

  A string naming the column in the condition slot that contains the
  group labels.

- n_components:

  The number of components for the PLS-DA model.

## Functions

- `plot_plsda(object = ProtData, group_variable = character)`: PLS-DA
  method for ProtData class.
