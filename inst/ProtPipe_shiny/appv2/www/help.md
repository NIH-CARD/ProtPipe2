# ProtPipe2 Help

ProtPipe2 is organized as a workflow:

1. Input
2. Quality Control
3. Pre-processing
4. Clustering / Dimensionality Reduction
5. Differential Intensity
6. Abundance Profiling

Use the tabs in order for a first pass through a dataset.

---

## Input

### Goal

Load a proteomics data matrix and, optionally, a sample metadata table.

### Accepted inputs

- CSV or TSV tables
- Excel files (`.xls`, `.xlsx`)
- SomaScan ADAT files
- Olink NPX files

### Intensity matrix

The main input table should contain:

- one row per protein or analyte
- one or more metadata columns describing each row
- one numeric column per sample

Example:

| ProteinID | Genes | Description | Sample_A | Sample_B | Sample_C | Sample_D |
|:---|:---|:---|---:|---:|---:|---:|
| P12345 | GENE1 | Protein 1 | 10234 | 9844 | 15322 | 14987 |
| P67890 | GENE2 | Protein 2 | 5531 | 6010 | 8122 | 7901 |
| Q11111 | GENE3 | Protein 3 | 0 | 214 | 6620 | 7012 |

### Sample metadata table

The metadata table is optional, but it is required for most grouped analyses.

Rules:

- it must contain a column named `SampleID`
- `SampleID` values must exactly match the sample column names in the intensity matrix
- additional columns can contain group labels, timepoints, batches, or continuous variables

Example:

| SampleID | Condition | Timepoint | Batch |
|:---|:---|:---|:---|
| Sample_A | Control | Day0 | Batch1 |
| Sample_B | Control | Day0 | Batch1 |
| Sample_C | Treated | Day28 | Batch2 |
| Sample_D | Treated | Day28 | Batch2 |

If no metadata table is provided, ProtPipe2 creates a default column called `base_condition` by removing a trailing sample index such as `_1` or `_2`.

---

## Quality Control

### Goal

Check sample quality, missingness, replicate consistency, and overall similarity before preprocessing.

### Main outputs

- **Protein Group Counts:** number of detected proteins per sample
- **Protein Intensity Plot:** distribution of intensities across samples
- **CV Plot:** variability of proteins within each condition
- **Sample Correlation Heatmap:** similarity between samples

### How to use this tab

- Use **Protein Group Counts** to find samples with unusually low coverage.
- Use **Protein Intensity Plot** to assess whether sample distributions are aligned.
- Use **CV Plot** to assess replicate consistency within groups.
- Use **Sample Correlation Heatmap** to confirm that replicates are more similar to each other than to other groups.

---

## Pre-processing

### Goal

Prepare the dataset for downstream analysis by filtering, transforming, normalizing, imputing, and correcting unwanted technical effects.

### Main actions

- remove proteins missing in too many samples
- remove outlier samples
- log2 transform intensities
- normalize by mean or median
- impute missing values
- apply batch correction

### How to use this tab

- Filter first if there are obvious low-quality proteins or samples.
- Log2 transform before most downstream analyses.
- Normalize if sample intensity distributions differ for technical reasons.
- Impute before PCA, UMAP, hierarchical clustering, or differential analysis if missing values remain.
- Use batch correction only when you have a known batch variable in the metadata table.

After preprocessing, return to **Quality Control** to confirm that the data look improved.

---

## Clustering / Dimensionality Reduction

### Goal

Visualize sample-level structure after preprocessing.

### Main outputs

- **PCA:** major axes of variation across samples
- **Hierarchical Clustering:** dendrogram of sample similarity
- **UMAP:** low-dimensional embedding that can separate non-linear structure

### How to interpret this tab

- Replicates should cluster together.
- Samples from distinct biological groups should separate when there is strong signal.
- If samples cluster by batch rather than condition, consider batch correction.

This tab typically works best after log transformation, normalization, and imputation.

---

## Differential Intensity

### Goal

Identify proteins associated with a categorical comparison or a continuous variable.

### Analysis modes

- **Two-group comparison:** compares one group against another
- **Continuous analysis:** identifies proteins correlated with a numeric variable

### Main outputs

- **Volcano plot:** effect size versus statistical significance
- **Results table:** statistics for each protein
- **GO pathway analysis:** biological interpretation of the differential signal

### How to interpret this tab

- Positive `logFC` values indicate higher abundance in the treatment group.
- Negative `logFC` values indicate lower abundance in the treatment group.
- Use adjusted p-values to identify statistically significant proteins.
- Use GO enrichment to summarize affected biological processes.

For pathway analysis, the selected gene symbol column must contain gene symbols that can be mapped to Entrez IDs.

---

## Abundance Profiling

### Goal

Inspect specific proteins or selected feature sets in more detail.

### Main outputs

- **Proteomics Heatmap:** visualizes a selected set of proteins across samples or conditions
- **Protein Barchart:** shows abundance of a single protein across samples or groups

### How to use this tab

- Use the heatmap for a targeted gene list or pathway-related set of proteins.
- Use the protein barchart to inspect a candidate biomarker or validate a differential result.

This tab is best used after differential analysis has identified proteins of interest.

---

## Data Privacy

Uploaded data is processed in-memory on Posit's cloud servers (AWS) and is not retained after your session ends. Data is transmitted over HTTPS and is not shared across user sessions. Do not upload data subject to data use agreements, IRB restrictions, or HIPAA regulations. For analyses involving sensitive or restricted data, use the [ProtPipe2 R package](https://github.com/NIH-CARD/ProtPipe2) locally.

---

## Common Issues

- Sample names in the intensity matrix and metadata table must match exactly.
- Clustering methods require missing values to be imputed first.
- Differential analysis requires at least two samples per group.
- Pathway analysis requires a gene symbol column that maps successfully to Entrez IDs.
