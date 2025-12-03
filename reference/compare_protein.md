# — 1. Required Packages —Make sure you have these packages installed for the function to workinstall.packages(c("dplyr", "ggplot2", "ggpubr", "rlang"))

— 1. Required Packages —Make sure you have these packages installed for
the function to workinstall.packages(c("dplyr", "ggplot2", "ggpubr",
"rlang"))

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

## — 2. Generic Function Definition —

\#' Generate a bar chart of protein intensity values \#' \#' Creates a
ggplot bar chart comparing the intensity of a single protein \#' either
across all samples or grouped by a condition. \#' \#' @param object The
`SummarizedExperiment` object. \#' @param prot A string specifying the
name of the protein to plot. \#' @param prot_meta_col A string naming
the column in the `rowData` slot to search for the protein. Defaults to
the first column. \#' @param condition (Optional) A string naming the
column in the `colData` slot to group samples by. \#' @return A ggplot
object. \#' @export setGeneric("compare_protein", def = function(object,
prot, prot_meta_col = NULL, condition = NULL)
standardGeneric("compare_protein") )

## — 3. S4 Method Implementation —
