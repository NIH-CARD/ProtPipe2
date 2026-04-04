# Run Over-Representation Analysis for Custom Gene Sets

A thin wrapper around
[`clusterProfiler::enricher()`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html)
for user-supplied gene sets represented as `TERM2GENE` and optional
`TERM2NAME` tables.

## Usage

``` r
enrich_terms(
  gene_id,
  all_gene_vector,
  term2gene,
  term2name = NULL,
  enrich_pvalue = 1
)
```

## Arguments

- gene_id:

  A character vector of significant gene IDs.

- all_gene_vector:

  A character vector containing the background universe.

- term2gene:

  A two-column data frame with term IDs in the first column and gene IDs
  in the second.

- term2name:

  An optional two-column data frame with term IDs and readable term
  names.

- enrich_pvalue:

  The p-value and q-value cutoff. Defaults to `1`.

## Value

An `enrichResult` object, or `NULL` if no terms are enriched.
