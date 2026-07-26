# ProtPipe2: proteomics analysis tool

Reproducible workflows for downstream proteomics analysis. All exported
functions accept and return `SummarizedExperiment` objects.

## Imported accessors

Package code calls the `SummarizedExperiment` accessors unqualified, so
they are imported here rather than repeated across the individual method
files. The same applies to the `stats`, `utils` and `dplyr` helpers used
throughout the package.

Note that `data.table` is deliberately *not* imported: doing so would
make the package "data.table aware" and change `[` semantics package
wide, so its functions are called with `data.table::` instead.

## See also

Useful links:

- <https://nih-card.github.io/ProtPipe2/>

## Author

**Maintainer**: Jacob Epstein <jacobepstein02@gmail.com>
([ORCID](https://orcid.org/0009-0003-7979-3532))

Authors:

- Jacob Epstein <jacobepstein02@gmail.com>
  ([ORCID](https://orcid.org/0009-0003-7979-3532))
