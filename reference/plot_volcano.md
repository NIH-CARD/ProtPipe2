# Plot a Volcano Plot for Differential Expression Results

This function generates a volcano plot based on log fold changes (logFC)
and adjusted p-values from differential expression analysis. It
highlights genes that pass the specified thresholds for logFC and FDR
(false discovery rate). Optionally, it can label genes of interest or
the top up/downregulated genes.

## Usage

``` r
plot_volcano(
  DT.original,
  label_col = NULL,
  lfc_threshold = 1,
  fdr_threshold = 0.01,
  labelgene = NULL
)
```

## Arguments

- DT.original:

  A data frame containing the differential expression results. It should
  have columns `logFC` (log fold change) and `adj.P.Val` (adjusted
  p-value).

- label_col:

  The column name (as a string) used for labeling genes in the plot.
  Default is `NULL`. If `NULL`, the function will select the first
  column of `DT.original` for labeling.

- lfc_threshold:

  A numeric value representing the threshold for log fold change
  (default is 1). Genes with logFC greater than or equal to this value
  are labeled as "UP", and genes with logFC less than or equal to the
  negative of this value are labeled as "DOWN".

- fdr_threshold:

  A numeric value representing the false discovery rate (FDR) threshold
  (default is 0.01). Genes with an adjusted p-value greater than or
  equal to this threshold will be labeled as "Others".

- labelgene:

  A character vector of gene names to be labeled in the plot (default is
  `NULL`). If provided, only these genes will be labeled in the plot.

## Value

A `ggplot2` object representing the volcano plot.
