# Protpipe Help & Tutorial

Welcome! This guide explains each panel of the Protpipe workflow and provides help on how to interpret the results.

---

## 1. Input Parameters

This is the first step of the analysis. You must upload your data and metadata here.

### Accepted File Formats
* **Text:** CSV (Comma-separated values) or TSV (Tab-separated values)
* **Excel:** `xls` or `xlsx`
* **SomaLogic:** `Somascan.adat` format

### Protein Intensity File

The main input is a dataframe containing your protein expression data.
This may contain a file from the SomaLogic or Olink proteomic platforms.
Otherwise it must be a tabular data file that must follow these rules:

- **Rows**: Each row should represent a single protein or analyte.
- **Columns**: The columns should contain both metadata about the
  proteins (e.g., Protein ID, Gene Name, UniProt Accession) and the
  quantitative values for each sample. The column headers for your
  samples will be used to link to the condition data.

### Example Protein Intensity File

Here is a small example of a correctly formatted protein data dataframe.
Columns `ProteinID` and `Gene.Name` are metadata, while `Sample_A`,
`Sample_B`, `Sample_C`, and `Sample_D` contain the measurements for each
sample.

| ProteinID | Gene.Name | Description | Sample_A | Sample_B | Sample_C | Sample_D |
|:---|:---|:---|---:|---:|---:|---:|
| P12345 | GENEA | Protein A Description | 1.2 | 1.5 | 5.5 | 5.9 |
| P67890 | GENEB | Protein B Description | 2.5 | 2.8 | 6.2 | 6.5 |
| Q54321 | GENEC | Protein C Description | 3.1 | 3.3 | 7.8 | 7.5 |
| Q09876 | GENED | Protein D Description | 4.0 | 4.2 | 8.1 | 8.3 |

### Sample Condition File

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

### Example Condition Format

This is an example of a correctly formatted condition dataframe that
corresponds to the protein data example above.

| SampleID | Condition | Timepoint | Batch |
|:---------|:----------|:----------|------:|
| Sample_A | Control   | 24h       |     1 |
| Sample_B | Control   | 24h       |     1 |
| Sample_C | Treated   | 24h       |     2 |
| Sample_D | Treated   | 24h       |     2 |

**Important**: A mismatch in sample names between the protein data
columns and the `SampleID` column in the condition file will result in
an error. Please ensure they are identical.

**Important**: If no Sample Condition File is provided, one will be 
generated internally. This will have one metadata column called "base_condition"
that contains the sample names with trailing "_integer" removed. For example if the
sample names are control_1, control_2, treatment_1, treatment_2, then the base_condition
column will contain control, control, treatment, treatment.

---

## 2. Quality Control (QC)

This panel helps you assess the quality of your data before analysis.

* **Coefficient of Variation (CV) per Protein:**
    * **What it is:** A measure of variability for each protein across all samples.
    * **How to interpret:** Lower CV is generally better. Proteins with very high CV may be unreliable.

* **Protein Intensity:**
    * **What it is:** Shows the distribution of protein intensities for each sample.
    * **How to interpret:** Before preprocessing, samples may have different distributions. After normalization (in the next step), you should return here and check that the distributions are properly aligned (e.g., all boxplot medians are at the same level).

* **Number of Non-Missing Proteins per Sample:**
    * **What it is:** A bar plot showing the total count of proteins detected in each sample.
    * **How to interpret:** Look for outliers. A sample with significantly fewer proteins than all others might be a failed run or a poor-quality sample that you may consider removing in the next step.

* **Spearman Correlation between Samples:**
    * **What it is:** A heatmap showing how similar your samples are to each other.
    * **How to interpret:** You want to see high correlation between your biological replicates. You should also see clear "blocks" of correlation corresponding to your experimental groups.

---

## 3. Preprocessing

This panel allows you to clean, transform, and correct your data.

* **Remove Samples / Remove Proteins:** 
    * **Remove outlier samples:** Removes samples (columns) if the number of non-missing proteins is outside n standard deviations from the mean
    * **Remove outlier proteins:** Removes proteins (rows) if they are present in less than n% of samples
* **Log2-Transform:**
    * **What it is:** Applies a `log_2(x+1)` transformation to your intensities.
    * **Why:** This is standard practice. It makes the data more symmetrical and prevents highly-abundant proteins from dominating the analysis.
* **Normalization:**
    * **What it is:** Corrects for technical variation between samples by leveling the global mean or median protein intensity (e.g., one sample was loaded with more total protein than another).
    * **How to interpret:** After normalizing, go back to the **QC** tab. The "Protein Intensity" boxplots should now be well-aligned.
* **Impute Missing Values:** Fill in missing data points (NAs). This is often required for clustering and differential.
    * **Fixed value:** Fill in missing values with n
    * **Minimum:** Fill in missing values with n(minimum value per protein)
    * **Left-shifted distribution:** Fits a normal distribution to impute missing values stochastically
* **Batch Correction:**
    * **Requirement:** You must have provided a metadata file.
    * **What it is:** Uses the `removeBatchEffect` function from the `limma` package to remove unwanted variation caused by processing samples in different batches.
    * **How to interpret:** After running this, check the **Clustering** tab. Samples should now not cluster by the batch variable.

---

## 4. Clustering / Dimensionality Reduction

This panel provides visual methods to see how your samples group together based on their protein expression.

* **Principal Component Analysis (PCA):**
    * **What it is:** A dimensionality-reduction plot. It condenses all protein information into two dimensions (PC1 and PC2) that capture the most variation.
    * **How to interpret:** Points (samples) that are closer together are more similar. You want to see your biological replicates clustering tightly. You should also see clear separation (distance) between your different experimental groups (e.g., "Control" vs. "Treated").

* **Hierarchical Clustering:**
    * **What it is:** A "tree" (dendrogram) that groups samples based on similarity.
    * **How to interpret:** Samples that join "lower" on the tree are more similar to each other. Like PCA, you want to see your replicates join first, and your main groups form distinct branches.

* **UMAP (Uniform Manifold Approximation and Projection):**
    * **What it is:** A more modern clustering method, often better at finding complex non-linear patterns than PCA.
    * **How to interpret:** Similar to PCA, look for distinct "islands" of samples that correspond to your biological groups.

---

## 5. Differential Expression (DE)

This panel identifies proteins that are significantly different between your experimental conditions. **Requires a metadata file.**

### Analysis Types
1.  **Categorical (2 Groups):**
    * **Use Case:** Comparing "Control" vs. "Treated".
    * **Method:** Uses the `limma` package to perform statistical testing.
2.  **Continuous Variable:**
    * **Use Case:** Finding proteins that correlate with a continuous variable, like "Age" or "Disease Score".
    * **Method:** Calculates Spearman's correlation.

### How to Interpret Results
* **Volcano Plot:**
    * **What it is:** A scatter plot of statistical significance (y-axis) vs. magnitude of change (x-axis).
    * **Interpretation:**
        * **Top Right (Significant & Up):** Proteins that are significantly up-regulated in your comparison.
        * **Top Left (Significant & Down):** Proteins that are significantly down-regulated.
        * **Middle/Bottom (Not Significant):** Proteins that are not statistically different.

* **DE Results Table:**
    * **`log2FoldChange` (log2FC):** The magnitude of change. A `log2FC` of 1 means a 2-fold increase. A `log2FC` of -2 means a 4-fold decrease.
    * **`pvalue`:** The raw statistical p-value.
    * **`padj` (or `FDR`):** The **adjusted p-value**. This is the most important value to use.
    * **Key Rule:** Proteins are typically considered **statistically significant** if their **`padj` is less than 0.05**.

* **Pathway Enrichment (ORA & GSEA):**
    * **What it is:** After finding *which* proteins are significant, this finds *which biological processes* (from GO and KEGG databases) are over-represented in your significant list.
    * **How to interpret:** This gives you biological context. For example, your results might show that "immune response pathways" or "metabolic pathways" are highly affected by your treatment.

---

## 6. Abundance Profiling

This panel allows you to "zoom in" on specific proteins of interest from your DE results or from prior knowledge.

* **Proteomics Heatmap:**
    * **What it is:** A grid of colors showing the expression of selected proteins (rows) across all your samples (columns).
    * **How to interpret:** This is a powerful way to see patterns. Look for "blocks" of color. For example, a red block might show a set of proteins that are all highly up-regulated in your "Treated" group, suggesting they are part of a coordinated biological response.

* **View Single Protein:**
    * **What it is:** Generates a barchart for any protein you select, either across all samples or between experimental groups. Pairwise comparison of means is done using the Wilcoxon Rank Sum Test.
    * **How to interpret:** This is a simple, direct way to visualize the expression of a single protein (e.g., "Protein_X") and see how it differs between your experimental groups.

### Saving Outputs
All plots generated in the app can be saved in **PDF format** for your records or for publication.
