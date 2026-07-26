# Perform Comprehensive GO and KEGG Pathway Enrichment Analysis

This function conducts a full suite of pathway enrichment analyses on a
given differential expression (DE) result set. It performs
Over-Representation Analysis (ORA) for significantly up- and
down-regulated genes, as well as Gene Set Enrichment Analysis (GSEA) on
the complete ranked list of genes. It analyzes both Gene Ontology (GO)
terms and KEGG pathways.

## Usage

``` r
enrich_pathways(
  DE,
  lfc_threshold = 1,
  fdr_threshold = 0.01,
  enrich_pvalue = 0.05,
  go_org = org.Hs.eg.db,
  kegg_org = "hsa",
  gene_col = "Genes",
  adj = TRUE,
  source = c("go", "custom"),
  go_ont = "BP",
  term2gene = NULL,
  term2name = NULL,
  run_ora = TRUE,
  run_gsea = TRUE,
  run_kegg = FALSE
)
```

## Arguments

- DE:

  A data frame containing differential expression results. Must include
  columns for gene identifiers, log fold change, and adjusted p-values.

- lfc_threshold:

  A numeric value for the absolute log2 fold change threshold to define
  significant genes for ORA. Defaults to `1`.

- fdr_threshold:

  A numeric value for the adjusted p-value (FDR) threshold to define
  significant genes for ORA. Defaults to `0.01`.

- enrich_pvalue:

  A numeric value for the p-value cutoff used to determine significant
  enrichment for pathways/terms. Defaults to `0.05`.

- go_org:

  An annotation database object (e.g., `org.Hs.eg.db` for human) used
  for GO analysis and mapping gene IDs.

- kegg_org:

  A character string specifying the KEGG organism code (e.g., `'hsa'`
  for Homo sapiens).

- gene_col:

  A character string indicating the name of the column in the `DE` data
  frame that contains the gene symbols/identifiers. Defaults to
  `"Genes"`.

- adj:

  Logical; when `TRUE` (the default) select significant genes using the
  adjusted p-value column, otherwise use the raw p-value.

- source:

  Pathway source. Use `"go"` for Gene Ontology or `"custom"` for
  user-supplied ontologies.

- go_ont:

  GO ontology to test when `source = "go"`. Defaults to `"BP"`.

- term2gene:

  Optional two-column data frame of term-to-gene mappings for custom
  ontologies. Gene IDs must match the IDs used for enrichment, currently
  Entrez IDs in the Shiny app workflow.

- term2name:

  Optional two-column data frame of term-to-name mappings for custom
  ontologies.

- run_ora:

  Logical; run over-representation analysis.

- run_gsea:

  Logical; run gene set enrichment analysis.

- run_kegg:

  Logical; also run KEGG analyses. Defaults to `FALSE`.

## Value

A list containing two named elements:

- `results`:

  A list of data frames with the detailed enrichment statistics for each
  analysis type (e.g., `go_up`, `kegg_down`, `gse_go`).

- `plots`:

  A list of `ggplot` objects for visualizing the enrichment results
  (e.g., dot plots, bar plots).

The function returns `NULL` if the initial mapping of gene identifiers
to Entrez IDs fails.

## Examples

``` r
# Create a sample DE results dataframe
de_results <- data.frame(
  Genes = c("TP53", "EGFR", "BRCA1", "TNF", "MMP9", "IL6", "VEGFA", "JUN"),
  logFC = c(2.5, -1.8, 2.1, 1.6, -2.2, 3.0, -1.7, 2.2),
  adj.P.Val = c(0.001, 0.005, 0.002, 0.009, 0.003, 0.0001, 0.008, 0.001)
)

# This example requires an internet connection and the org.Hs.eg.db package.
if (FALSE) { # \dontrun{
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  enrichment_output <- enrich_pathways(
    DE = de_results,
    lfc_threshold = 1.5,
    fdr_threshold = 0.01,
    go_org = org.Hs.eg.db,
    kegg_org = 'hsa'
  )

  # Check if the analysis produced results
  if (!is.null(enrichment_output)) {
    # View the head of the GO results for upregulated genes
    print(head(enrichment_output$results$go_up))

    # Display one of the generated plots
    if (!is.null(enrichment_output$plots$go_up_dotplot)) {
      print(enrichment_output$plots$go_up_dotplot)
    }
  }
}
} # }
```
