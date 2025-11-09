# Add Entrez Gene IDs to a Data Frame

Maps gene symbols from a specified column to Entrez IDs using
clusterProfiler::bitr and joins them back to the original data frame. It
also cleans gene symbols that may contain multiple entries separated by
semicolons.

## Usage

``` r
add_entrez(DE, org = org.Hs.eg.db::org.Hs.eg.db, gene_col = "Genes")
```

## Arguments

- DE:

  A data frame, such as one containing differential expression results.

- org:

  An organism annotation database object (e.g., org.Hs.eg.db).

- gene_col:

  A character string specifying the name of the column in DE that
  contains the gene symbols.

## Value

A data frame with an appended ENTREZID column containing the mapped IDs.
Returns NULL if no genes can be successfully mapped.
