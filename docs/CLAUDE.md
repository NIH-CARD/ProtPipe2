# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## What This Is

ProtPipe2 is an R package (`ProtPipe`) for downstream proteomics
analysis, paired with a Shiny web app. It supports data import into
`SummarizedExperiment`, quality control, preprocessing, dimensionality
reduction, clustering, differential intensity analysis, and pathway
analysis.

## Common Commands

Run from an R console in the project root.

``` r

# Load package for development (run after any change to R/)
devtools::load_all()

# Run all tests
devtools::test()

# Run a single test file
testthat::test_file("tests/testthat/test-preprocessing_methods.R")

# Regenerate man/ documentation from roxygen comments
devtools::document()

# Full R CMD check
devtools::check()

# Install the package
devtools::install()

# Launch the Shiny app (after installing)
ProtPipe::run_protpipe_shiny()
```

## Architecture

### Data Model

`SummarizedExperiment` is the only supported data structure. All
exported package functions accept and return `SummarizedExperiment`
objects. The `ProtData` class was a previous implementation — do not
reintroduce it.

- `assay(se)` — the intensity matrix (proteins × samples)
- `rowData(se)` — protein metadata (gene names, protein IDs, etc.)
- `colData(se)` — sample metadata (conditions, batch, covariates)
- `metadata(se)$processing_log` — ordered list of applied preprocessing
  steps, each a named list with `name`, `parameters`, and `details`

The `has_step(se, "step_name")` helper checks whether a preprocessing
step has been applied, which functions like `do_limma_binary` use to
conditionally apply `log2_transform`.

### R/ Source Files

| File | Purpose |
|----|----|
| `AllGenerics.R` | All S4 generic definitions (the single source of truth for function signatures) |
| `0-prot_data.R` | [`create_se()`](https://nih-card.github.io/ProtPipe2/docs/reference/create_se.md) constructor, [`has_step()`](https://nih-card.github.io/ProtPipe2/docs/reference/has_step.md), [`detect_intensity_cols()`](https://nih-card.github.io/ProtPipe2/docs/reference/detect_intensity_cols.md), internal helpers |
| `preprocessing_methods.R` | S4 method implementations: filtering, normalization, imputation, batch correction; also [`add_processing_step()`](https://nih-card.github.io/ProtPipe2/docs/reference/add_processing_step.md) and [`generate_preprocessing_report()`](https://nih-card.github.io/ProtPipe2/docs/reference/generate_preprocessing_report.md) |
| `quality_control.R` | QC plots and metrics (protein counts, intensity distributions, CVs, sample correlation) |
| `clustering.R` | PCA, UMAP, hierarchical clustering |
| `differential_expression.R` | limma, t-test, Spearman correlation; volcano plots; pathway analysis (GO/KEGG/custom ontologies) |
| `heatmap.R` | [`plot_proteomics_heatmap()`](https://nih-card.github.io/ProtPipe2/docs/reference/plot_proteomics_heatmap.md) |
| `protein_comparison.R` | [`compare_protein()`](https://nih-card.github.io/ProtPipe2/docs/reference/compare_protein.md) for abundance profiling |
| `olink.R` | [`create_se_from_olink()`](https://nih-card.github.io/ProtPipe2/docs/reference/create_se_from_olink.md) for Olink NPX data import |
| `soma.R` | [`create_se_from_soma()`](https://nih-card.github.io/ProtPipe2/docs/reference/create_se_from_soma.md) for SomaScan ADAT import |
| `run_shiny.R` | [`run_protpipe_shiny()`](https://nih-card.github.io/ProtPipe2/docs/reference/run_protpipe_shiny.md) launcher |
| `dependencies.R` | `protpipe_require_packages()` helper for optional Suggests dependencies |

### Generics Pattern

All exported generics are declared in `AllGenerics.R` using
`setGeneric()`. Method implementations live in the domain-specific files
using `setMethod()`. When adding a new exported function, add its
generic to `AllGenerics.R` first.

### Shiny App (`inst/ProtPipe_shiny/`)

The app is a thin wrapper around package functions — analysis logic
belongs in the package, not the Shiny server. The workflow tabs match:
Input → Quality Control → Pre-Processing → Clustering / Dimensionality
Reduction → Differential Intensity → Abundance Profiling → Help

Help content is in `inst/ProtPipe_shiny/www/help.md`; if it changes,
regenerate `help.html`.

### Data Import Paths

- **Generic tabular data**: `create_se(data, sample_metadata)` —
  autodetects numeric columns as intensities; `SampleID` column required
  in metadata
- **Olink**: `create_se_from_olink(npx, condition)` — wraps
  [`OlinkAnalyze::read_npx()`](https://rdrr.io/pkg/OlinkAnalyze/man/read_npx.html)
  output
- **SomaScan**: `create_se_from_soma(adat, condition)` — wraps
  `SomaDataIO` ADAT objects

### Testing

Tests live in `tests/testthat/`. `helper-fixtures.R` provides
`load_basic_se()` and related helpers that read from
`EXAMPLES/basic_example_data/`. Use relative paths via `package_path()`
— never [`setwd()`](https://rdrr.io/r/base/getwd.html) or absolute paths
in tests.

### Optional Dependencies

Heavy Bioconductor packages (org.Hs.eg.db, OlinkAnalyze, SomaDataIO,
DelayedMatrixStats) are in `Suggests`. The
`protpipe_require_packages("pkg", feature = "fn()")` helper in
`dependencies.R` provides a user-friendly error when they are missing.
