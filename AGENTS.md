# AGENTS.md

## Project Scope

This repository contains:

- an R package for downstream proteomics analysis in `R/`
- a Shiny app interface in `inst/ProtPipe_shiny/`

The package is the primary source of analytical logic. The Shiny app should call package functions rather than duplicating analysis code.

## Data Model

- Use `SummarizedExperiment` as the only supported core data structure.
- Do not reintroduce `ProtData` or add new code paths that depend on it.
- New package functions should accept and return `SummarizedExperiment` objects unless there is a strong reason not to.

## Package Conventions

- Document all new exported functions with roxygen comments.
- Prefer updating existing package functions over adding Shiny-only logic.
- Keep function interfaces simple and consistent with the current package style.
- When editing docs or examples, remove stale references rather than preserving backwards-looking commentary.

## Shiny Conventions

- Keep the Shiny app as a thin wrapper around package functionality.
- Put parsing, validation, and analysis logic in the package where practical.
- Avoid adding display-only helper columns or other superfluous UI embellishments that were not explicitly requested.
- Preserve the workflow structure:
  - Input
  - Quality Control
  - Pre-Processing
  - Clustering / Dimensionality Reduction
  - Differential Intensity
  - Abundance Profiling
  - Help
- Keep labels user-facing and plain English. Do not expose internal function names in the UI.

## Testing

- Add or update tests for package-facing behavior when changing analytical code.
- Prefer bundled example data in `EXAMPLES/` over absolute paths or local machine files.
- Use root-aware test paths. Do not rely on `setwd()` in tests.
- Smoke tests are acceptable for plotting and app-adjacent helpers.

## Documentation Assets

- Shiny help content lives in `inst/ProtPipe_shiny/www/help.md`.
- If `help.md` changes, regenerate `inst/ProtPipe_shiny/www/help.html`.

## Benchmarks

- Benchmark scripts live in `benchmarks/`.
- Keep benchmark utilities out of `inst/`.
- Prefer OS-level memory measurement over in-R peak memory estimates when reporting compute usage.

## Cleanup Expectations

- Avoid leaving dead legacy files, stale examples, or references to removed APIs.
- If a feature is removed, update tests, roxygen comments, and UI text to match.
