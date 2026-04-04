# Perform Gene Ontology (GO) Gene Set Enrichment Analysis (GSEA)

A wrapper around clusterProfiler::gseGO to perform GSEA for Gene
Ontology terms on a ranked list of genes.

## Usage

``` r
gse_go(
  gene_list,
  enrich_pvalue = 1,
  org = org.Hs.eg.db::org.Hs.eg.db,
  ont = "BP"
)
```

## Arguments

- gene_list:

  A named numeric vector of genes. Names should be Entrez IDs and values
  should be the ranking metric (e.g., log2 fold change).

- enrich_pvalue:

  The p-value cutoff for enrichment. Defaults to 1 (no cutoff).

- org:

  An organism annotation database (e.g., org.Hs.eg.db).

- ont:

  GO ontology to test. One of `"BP"`, `"MF"`, `"CC"`, or `"ALL"`.

## Value

A gseResult object if the analysis is successful and finds enriched
terms, otherwise NULL.
