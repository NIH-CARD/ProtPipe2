# Run GSEA for Custom Gene Sets

A thin wrapper around
[`clusterProfiler::GSEA()`](https://rdrr.io/pkg/clusterProfiler/man/GSEA.html)
for user-supplied gene sets represented as `TERM2GENE` and optional
`TERM2NAME` tables.

## Usage

``` r
gse_terms(gene_list, term2gene, term2name = NULL, enrich_pvalue = 1)
```

## Arguments

- gene_list:

  A named numeric vector of ranked genes.

- term2gene:

  A two-column data frame with term IDs in the first column and gene IDs
  in the second.

- term2name:

  An optional two-column data frame with term IDs and readable term
  names.

- enrich_pvalue:

  The p-value cutoff. Defaults to `1`.

## Value

A `gseaResult` object, or `NULL` if no terms are enriched.
