Input Data Format
================

This page details the required format for the input data used in this R
package. The package requires one primary data file (for protein/analyte
values) and can optionally accept a second file containing experimental
conditions or metadata.

## 1. Protein Data File

The main input is a dataframe containing your protein expression data.
The structure must follow these rules:

- **Rows**: Each row should represent a single protein or analyte.
- **Columns**: The columns should contain both metadata about the
  proteins (e.g., Protein ID, Gene Name, UniProt Accession) and the
  quantitative values for each sample. The column headers for your
  samples will be used to link to the condition data.

### Example Protein Data Table

Here is a small example of a correctly formatted protein data dataframe.
Columns `ProteinID` and `Gene.Name` are metadata, while `Sample_A`,
`Sample_B`, `Sample_C`, and `Sample_D` contain the measurements for each
sample.

| ProteinID | Gene.Name | Description           | Sample_A | Sample_B | Sample_C | Sample_D |
|:----------|:----------|:----------------------|---------:|---------:|---------:|---------:|
| P12345    | GENEA     | Protein A Description |      1.2 |      1.5 |      5.5 |      5.9 |
| P67890    | GENEB     | Protein B Description |      2.5 |      2.8 |      6.2 |      6.5 |
| Q54321    | GENEC     | Protein C Description |      3.1 |      3.3 |      7.8 |      7.5 |
| Q09876    | GENED     | Protein D Description |      4.0 |      4.2 |      8.1 |      8.3 |

Example Protein Data Structure

------------------------------------------------------------------------

## 2. Condition Metadata File (Optional)

You can provide an optional second dataframe that describes the
experimental conditions for each sample. This file is highly recommended
for downstream analysis.

- **Crucial Requirement**: This dataframe **MUST** contain a column
  named `SampleID`.
- The values in the `SampleID` column **MUST** exactly match the column
  names of the samples in your protein data file (e.g., `Sample_A`,
  `Sample_B`, etc.).
- Other columns can contain any metadata you wish to associate with your
  samples, such as treatment group, time point, batch, etc.

### Example Condition Data Table

This is an example of a correctly formatted condition dataframe that
corresponds to the protein data example above.

| SampleID | Condition | Timepoint | Batch |
|:---------|:----------|:----------|------:|
| Sample_A | Control   | 24h       |     1 |
| Sample_B | Control   | 24h       |     1 |
| Sample_C | Treated   | 24h       |     2 |
| Sample_D | Treated   | 24h       |     2 |

Example Condition/Metadata Structure

**Important**: A mismatch in sample names between the protein data
columns and the `SampleID` column in the condition file will result in
an error. Please ensure they are identical.
