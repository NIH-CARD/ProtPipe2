# Creates a ggplot bar chart comparing the intensity of a single protein either across all samples or grouped by a condition.

Creates a ggplot bar chart comparing the intensity of a single protein
either across all samples or grouped by a condition.

## Usage

``` r
compare_protein(
  object,
  prot,
  prot_meta_col = NULL,
  condition = NULL,
  selected_groups = NULL
)

# S4 method for class 'SummarizedExperiment'
compare_protein(
  object,
  prot,
  prot_meta_col = NULL,
  condition = NULL,
  selected_groups = NULL
)
```

## Arguments

- object:

  The `SummarizedExperiment` object.

- prot:

  A string specifying the name of the protein to plot.

- prot_meta_col:

  A string naming the column in the `rowData` slot to search for the
  protein. Defaults to the first column.

- condition:

  (Optional) A string naming the column in the `colData` slot to group
  samples by.

- selected_groups:

  (Optional) A character vector of group names to filter by. Only used
  if 'condition' is provided.

## Value

A ggplot object.

## Functions

- `compare_protein(SummarizedExperiment)`: Method for
  SummarizedExperiment objects.
