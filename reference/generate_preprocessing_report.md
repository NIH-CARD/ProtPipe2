# Generate Markdown Report of Preprocessing Steps

Takes a SummarizedExperiment object with a populated `processing_log` in
its metadata and writes a Markdown file summarizing the preprocessing
steps in order.

## Usage

``` r
generate_preprocessing_report(object, output_file = "preprocessing_report.md")
```

## Arguments

- object:

  A SummarizedExperiment object.

- output_file:

  A character string giving the name of the output Markdown file.

## Value

Invisibly returns the file path to the Markdown report.
