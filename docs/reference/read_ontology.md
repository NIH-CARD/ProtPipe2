# Read a Custom Ontology File

Reads a user-supplied ontology into the `TERM2GENE` and optional
`TERM2NAME` structures required by
[`clusterProfiler::enricher()`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html)
and
[`clusterProfiler::GSEA()`](https://rdrr.io/pkg/clusterProfiler/man/GSEA.html).

Supported formats are:

- GMT files, where each line is `term`, `name`, then one or more gene
  IDs

- Delimited text files (`csv` or `tsv`) with either:

  - 2 columns: `term`, `gene`

  - 3 columns: `term`, `name`, `gene`

Gene identifiers in uploaded ontologies must match the identifiers used
for enrichment. In the ProtPipe Shiny app, this currently means Entrez
gene IDs.

## Usage

``` r
read_ontology(file)
```

## Arguments

- file:

  Path to a `.gmt`, `.csv`, or `.tsv` ontology file.

## Value

A named list with `term2gene` and `term2name`.
