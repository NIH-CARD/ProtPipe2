# Perform KEGG Over-Representation Analysis (ORA)

A wrapper around clusterProfiler::enrichKEGG to perform ORA to find
enriched KEGG pathways from a list of significant genes.

## Usage

``` r
enrich_kegg(
  gene_id,
  all_gene_vector,
  enrich_pvalue = 1,
  org = org.Hs.eg.db::org.Hs.eg.db,
  organism = "hsa"
)
```

## Arguments

- gene_id:

  A character vector of significant Entrez gene IDs to test for
  enrichment.

- all_gene_vector:

  A character vector of all background (universe) Entrez gene IDs
  against which to test.

- enrich_pvalue:

  The p-value and q-value cutoff for enrichment. Defaults to 1.

- org:

  An organism annotation database (e.g., org.Hs.eg.db) used to map
  Entrez IDs to gene symbols.

- organism:

  A character string for the KEGG organism code (e.g., 'hsa' for human).

## Value

An enrichResult object if the analysis is successful and finds enriched
pathways, otherwise NULL.
