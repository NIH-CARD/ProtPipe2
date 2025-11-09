# Perform Gene Ontology (GO) Over-Representation Analysis (ORA)

A wrapper around clusterProfiler::enrichGO to perform ORA to find
enriched GO terms (Biological Process, Molecular Function, Cellular
Component).

## Usage

``` r
enrich_go(
  gene_id,
  all_gene_vector,
  enrich_pvalue = 1,
  org = org.Hs.eg.db::org.Hs.eg.db
)
```

## Arguments

- gene_id:

  A character vector of significant Entrez gene IDs to test.

- all_gene_vector:

  A character vector of all background (universe) Entrez gene IDs
  against which to test.

- enrich_pvalue:

  The p-value and q-value cutoff for enrichment. Defaults to 1.

- org:

  An organism annotation database (e.g., org.Hs.eg.db).

## Value

An enrichResult object if the analysis is successful and finds enriched
terms, otherwise NULL.
