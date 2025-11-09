# Perform KEGG Gene Set Enrichment Analysis (GSEA)

A wrapper around clusterProfiler::gseKEGG to perform GSEA on a ranked
list of genes and convert IDs to readable gene symbols.

## Usage

``` r
gse_kegg(
  gene_list,
  enrich_pvalue = 1,
  org = org.Hs.eg.db::org.Hs.eg.db,
  organism = "hsa"
)
```

## Arguments

- gene_list:

  A named numeric vector of genes. Names should be Entrez IDs and values
  should be the ranking metric (e.g., log2 fold change).

- enrich_pvalue:

  The p-value cutoff for enrichment. Defaults to 1 (no cutoff).

- org:

  An organism annotation database (e.g., org.Hs.eg.db) used to map
  Entrez IDs to gene symbols.

- organism:

  A character string for the KEGG organism code (e.g., 'hsa' for human).

## Value

A gseKEGGResult object if the analysis is successful and finds enriched
pathways, otherwise NULL.
