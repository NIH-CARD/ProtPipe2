## Quality Control ############################################################

#' Get Protein Group Counts Per Sample
#'
#' @description
#' Calculates the number of non-missing protein groups (or features) for each
#' sample in a `SummarizedExperiment` object.
#'
#' @param object A \code{SummarizedExperiment} object.
#'
#' @return A data frame with two columns: `Sample` (the sample names) and `Protein_Groups`
#'   (the count of non-NA proteins for that sample).
#'
#' @export
#'
setGeneric("get_pg_counts", function(object) standardGeneric("get_pg_counts"))

#' Plot Protein Group Counts
#'
#' @description
#' Generates a bar plot showing the number of identified protein groups per sample.
#' The plot can either show counts for each individual sample, or show the
#' summarized mean and standard deviation for samples grouped by a condition.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param condition Optional. A character string specifying a column name in the
#'   `colData` slot. If provided, samples will be grouped by this variable and
#'   the plot will show the mean and standard deviation of counts per group.
#'   If `NULL` (the default), each sample is plotted individually.
#'
#' @return A `ggplot` object.
#'
#' @export
#'
setGeneric("plot_pg_counts", function(object, condition = NULL) standardGeneric("plot_pg_counts"))

#' Plot Boxplots of Sample Intensity Distributions
#'
#' @description
#' Generates boxplots to visualize and compare the distribution of protein
#' intensities for each sample. This is a common quality control plot to check
#' for normalization issues or identify potential outlier samples.
#'
#' @param object A \code{SummarizedExperiment} object.
#'
#' @return A `ggplot` object showing a boxplot for each sample's
#'   log10-transformed intensity distribution.
#'
#' @export
#'
#' @examples
#' # Create sample data
#' raw_data <- data.frame(
#'   Gene = c("GENEA", "GENEB", "GENEC", "GENED"),
#'   SampleA = c(100, 200, 150, 120),
#'   SampleB = c(110, 210, 160, 130), # Slightly higher than A
#'   SampleC = c(250, 500, 400, 300)   # Higher median and spread
#' )
#' se <- create_se(raw_data)
#'
#' # Generate the boxplot of intensities
#' p <- plot_pg_intensities(se)
#' if (interactive()) {
#'   print(p)
#' }
#'
setGeneric("plot_pg_intensities", function(object) standardGeneric("plot_pg_intensities"))

#' Calculate Coefficient of Variation (CV) for Protein Groups
#'
#' @description
#' Calculates the coefficient of variation (CV) for each protein, grouped by an
#' experimental condition. The CV is a standardized measure of dispersion,
#' calculated as the standard deviation divided by the mean, and is often used
#' to assess the reproducibility of replicates within a condition group.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param condition A character string specifying the column name in the
#'   `colData` slot which defines the replicate groups.
#' @param min_samples The minimum number of replicates in a group required to
#'   calculate CVs. Groups with fewer samples will be skipped. Defaults to 2.
#'
#' @return A data frame in long format with columns for `Protein`, `CV`
#'   (in percent), and `Condition`.
#'
#' @export
#'
setGeneric("get_CVs", function(object, condition, min_samples = 2) standardGeneric("get_CVs"))

#' Plot Coefficient of Variation (CV) Distributions
#'
#' @description
#' Generates a violin or jitter plot to visualize the distribution of protein
#' CVs for each experimental condition. This plot is useful for comparing the
#' reproducibility and technical variance across different sample groups.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param condition A character string specifying the column name in the
#'   `colData` slot which defines the replicate groups.
#' @param plot_type The type of plot to generate. Must be either `"violin"`
#'   (the default, which also includes a narrow boxplot) or `"jitter"`.
#'
#' @return A `ggplot` object.
#'
#' @export
#'
setGeneric("plot_CVs",
           function(object, condition, plot_type = "violin") {
             standardGeneric("plot_CVs")
           }
)

#' Calculate Pairwise Sample Correlations
#'
#' @description
#' Computes a pairwise correlation matrix for all samples based on their protein
#' intensity profiles. The result is returned in a long "tidy" format data table.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param method The correlation method to use. Can be one of `"spearman"` (the
#'   default), `"pearson"`, or `"kendall"`.
#' @param num_features Optional. The number of most variable proteins to use for
#'   the correlation calculation. If `NULL` (the default), all proteins are
#'   used.
#'
#' @return A `data.table` with three columns: `SampleA`, `SampleB`, and a third
#'   column named for the method used (e.g., `Spearman`), containing the
#'   pairwise correlation coefficients.
#'
#' @export
#'
setGeneric(
  "get_sample_correlation",
  function(object, method = "spearman", num_features = NULL) {
    standardGeneric("get_sample_correlation")
  }
)

#' Plot a Sample Correlation Heatmap
#'
#' @description
#' Creates a heatmap of pairwise sample correlations, calculated via the Spearman
#' method. The function includes a "smart sorting" feature to correctly order
#' axes based on numeric parts of labels (e.g., "Day1", "Day2", "Day10").
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param order_by Optional. A character string specifying a column in the
#'   condition slot by which to order the samples on the heatmap axes. If the
#'   column contains numeric parts, a natural sort is applied.
#' @param label_by Optional. A character string specifying a column in the
#'   condition slot to use for relabeling the heatmap axes.
#' @param num_features Optional. The number of most variable proteins to use for
#'   the correlation calculation. If `NULL` (the default), all proteins are
#'   used.
#'
#' @return A `ggplot` object representing the heatmap.
#'
#' @export
#'
setGeneric("plot_correlation_heatmap",
           function(object, order_by = NULL, label_by = NULL, num_features = NULL) {
             standardGeneric("plot_correlation_heatmap")
           }
)

## preprocessing ############################################################

#' Filter Assay Based on Limit of Detection (LOD)
#'
#' Filters the assay data of a SummarizedExperiment by comparing intensities
#' against a specific Limit of Detection (LOD). The LOD values can be sourced
#' either from a metadata column (rowData) or a specific reference sample
#' (e.g., a Buffer column in the assay).
#'
#' @param se A \code{SummarizedExperiment} object.
#' @param lod_col A character string indicating where to find the LOD values.
#'   The function first searches \code{rowData(se)}, then searches the
#'   column names of the assay (e.g., a specific buffer sample).
#'   Defaults to "Buffer".
#'
#' @return A \code{SummarizedExperiment} object where values below the
#'   defined LOD are replaced with \code{NA}.
#'
#' @export
setGeneric("lod_filter", function(se, lod_col = "Buffer") {
  standardGeneric("lod_filter")
})

#' Apply Limit of Detection Threshold
#'
#' Filters assay data by converting values below a specified Limit of Detection
#' (LOD) to \code{NA}.
#'
#' @param object A \code{SummarizedExperiment} object containing a single
#'   proteomics assay.
#' @param lod A single numeric value representing the Limit of Detection.
#'   Intensities below this value will be replaced with \code{NA}.
#'
#' @return The original \code{object} with the assay matrix updated to contain
#'   \code{NA} values where intensities were below the \code{lod}.
#'
#' @export
#'
#' @examples
#' # Create a mock SummarizedExperiment
#' counts <- matrix(c(10, 5, 2, 8, 1, 15), nrow = 3, ncol = 2)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#'
#' # Apply LOD of 6
#' se_filtered <- applyLOD(se, lod = 6)
#' assay(se_filtered)
setGeneric("apply_min_intenisty", function(object, lod) {
  standardGeneric("apply_min_intenisty")
})

#' @title Filter Proteins by Percentage of Valid Values
#'
#' @description
#' A generic function to filter proteins (rows) from a proteomics data object
#' based on the percentage of samples where they have valid observations.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param percent The minimum percentage (from 0 to 100) of samples in which a
#'   protein must have a valid (non-NA) value to be retained.
#'
#' @return A modified object of the same class with proteins filtered.
#'
#' @export
setGeneric("filter_proteins_by_percent",
           def = function(object, percent) {
             standardGeneric("filter_proteins_by_percent")
           }
)

#' Removes duplicate analytes
#'
#' Given a prot_meta column name, it will remove rows so that each value in the column is unique.
#' For each row the analyte with 1) the lowest missing values and 2) the greatest median intensity will be kept.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param col a column name of the rowData slot
#'
#' @return A \code{SummarizedExperiment} object.
#' @export
#'
setGeneric("filter_unique_proteins", function(object, col = NULL) standardGeneric("filter_unique_proteins"))

#' Filter Outlier Samples Based on Protein Counts
#'
#' @description
#' This function identifies and removes entire samples that are considered outliers
#' based on the total number of non-missing proteins identified within them.
#'
#' The method calculates the mean and standard deviation of protein counts across
#' all samples. A sample is flagged as an outlier if its total protein count
#' falls outside the range defined by a specified number of standard deviations
#' (`sds`) from the mean.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param sds The number of standard deviations from the mean protein count to use
#'   as the outlier threshold. Defaults to 3.
#'
#' @return A \code{SummarizedExperiment} object with the outlier samples removed.
#'
#' @export
#'
setGeneric("filter_outlier_samples", function(object, sds = 3) standardGeneric("filter_outlier_samples"))

#' Retain proteins present in a specified group of a \code{SummarizedExperiment} object
#'
#' Filters the object to retain only proteins that are present (non-NA)
#' in at least one sample within each unique group of the specified condition.
#'
#' @param object A \code{SummarizedExperiment} object
#' @param condition_name A character string specifying the column name in the
#' `colData` slot to group samples by.
#'
#' @return A new, modified \code{SummarizedExperiment} object containing only the overlapping proteins.
setGeneric("filter_overlap",
           def = function(object, condition_name) {
             standardGeneric("filter_overlap")
           }
)

#' Z-Score Normalization for Proteins Across Samples
#'
#' @description
#' This method performs a Z-score transformation on a protein-wise (row-wise)
#' basis. For each protein, it calculates the mean and standard deviation of its
#' abundance across all samples and then scales the values accordingly.
#'
#' This is a standard transformation for visualizing expression patterns in a
#' heatmap, as it highlights the relative change of each protein across samples,
#' independent of its absolute abundance.
#'
#' @param object A \code{SummarizedExperiment} object
#'
#' @return object A \code{SummarizedExperiment} object where the abundance values in the `assay` slot
#'   have been replaced by their row-wise Z-scores.
#'
#' @export
#'
#' @seealso [base::scale()]
#'
#' @examples
#' # Create sample data with proteins as rows and samples as columns
#' df <- data.frame(
#'   Protein = c("Protein1", "Protein2", "Protein3"),
#'   SampleA = c(100, 250, 50),
#'   SampleB = c(120, 200, 100),
#'   SampleC = c(110, 225, 75)
#' )
#'
#' sample_info <- data.frame(
#'   SampleID = c("SampleA", "SampleB", "SampleC"),
#'   group = c("Control", "Treatment", "Control")
#' )
#'
#' se <- create_se(df, sample_metadata = sample_info)
#'
#' # Check the means of each protein (row) before scaling
#' rowMeans(assay(se))
#'
#' # Apply the scaling method
#' scaled_se <- z_score(se)
#'
#' # The new data has row means near zero and row standard deviations of one
#' print(assay(scaled_se))
#' cat("Row means after scaling:\n")
#' print(rowMeans(assay(scaled_se)))
#' cat("\nRow standard deviations after scaling:\n")
#' print(apply(assay(scaled_se), 1, sd))
#'
setGeneric("z_score", function(object) standardGeneric("z_score"))

#' Median Normalization of Proteomics Data
#'
#' @description
#' Performs median normalization on the quantitative data within a
#' `SummarizedExperiment`
#' object. This method corrects for systematic, sample-specific biases (e.g.,
#' differences in sample loading or instrument sensitivity) to make the samples
#' more comparable.
#'
#' The function operates by first calculating the median intensity for each
#' sample (column). It then calculates the median of these medians (the "global
#' median"). Finally, it scales the intensities in each sample by a
#' multiplicative factor so that every sample has the same median abundance,
#' equal to the global median.
#'
#' @details
#' This type of multiplicative adjustment is typically appropriate for raw,
#' non-log-transformed intensity data.
#'
#' @param object A \code{SummarizedExperiment} object
#'
#' @return A \code{SummarizedExperiment} object with the abundance data in the `assay` slot
#'   normalized by the median-centering method.
#'
#' @export
setGeneric("median_normalize", function(object) standardGeneric("median_normalize"))

#' Mean Normalization of Proteomics Data
#'
#' @description
#' Corrects for sample-specific biases by scaling each sample (column) to have
#' the same mean abundance.
#'
#' This method adjusts each sample's intensities by a multiplicative factor,
#' aligning the mean of every sample to a global mean calculated from the
#' entire dataset.
#'
#' @details
#' This multiplicative adjustment is typically used for raw, non-log-transformed
#' intensity data.
#'
#' @param object A \code{SummarizedExperiment} object
#'
#' @return A \code{SummarizedExperiment} object with its data normalized by the mean-centering
#'   method.
#'
#' @export
setGeneric("mean_normalize", function(object) standardGeneric("mean_normalize"))

#' Performs a log2 transform of protein intensity values
#'
#' @param object A \code{SummarizedExperiment} object
#'
#' @return A \code{SummarizedExperiment} object
#' @export
#'
setGeneric("log2_transform", function(object) standardGeneric("log2_transform"))

#' Impute Missing Values with a Constant
#'
#' @description
#' Replaces all missing values (`NA` and `NaN`) in the numeric columns of the
#' data with a single, user-specified constant.
#'
#' @param object A \code{SummarizedExperiment} object containing data with missing values.
#' @param value The numeric constant to use for replacing `NA` and `NaN` values.
#'
#' @return A \code{SummarizedExperiment} object with missing values imputed.
#'
#' @export
#' @examples
#' # Create a sample data frame with metadata and a missing numeric value
#' raw_data <- data.frame(
#'   Protein.ID = c("P02768", "P01023", "P60709"),
#'   Gene.Name = c("ALB", "A2M", "ACTB"),
#'   Sample_A = c(1.2e6, 2.3e6, NA),
#'   Sample_B = c(1.4e6, 2.6e6, 4.8e6)
#' )
#'
#' # Create a SummarizedExperiment object
#' se <- create_se(raw_data)
#'
#' # Impute the NA with 0
#' imputed_obj <- impute(se, value = 0)
#'
#' # View the imputed assay data
#' print(assay(imputed_obj))
#'
setGeneric("impute", function(object, value) standardGeneric("impute"))

#' Impute Missing Values with the Row Minimum
#'
#' @description
#' Performs row-wise minimum imputation. For each protein (row), it replaces
#' missing values (`NA`, `NaN`) with the minimum observed value found in that
#' same row.
#'
#' The imputation value can be scaled by a multiplicative factor `alpha`.
#'
#' @param object A \code{SummarizedExperiment} object containing data with missing values.
#' @param alpha A numeric scaling factor to multiply the row minimum by before
#'   imputation. Defaults to 1 (no scaling).
#'
#' @return A \code{SummarizedExperiment} object with missing values imputed on a per-protein basis.
#'
#' @export
#' @examples
#' # Create data with different minimums and NAs in each row
#' raw_data <- data.frame(
#'   Gene = c("GENEA", "GENEB"),
#'   SampleA = c(100, 500),
#'   SampleB = c(200, 600),
#'   SampleC = c(NA, NA)
#' )
#'
#' se <- create_se(raw_data)
#' cat("Original Data:\n")
#' print(assay(se))
#'
#' # Impute using the row minimum (alpha = 1)
#' # Row 1's NA becomes 100; Row 2's NA becomes 500.
#' imputed_obj <- impute_min(se)
#' cat("\nImputed with alpha = 1:\n")
#' print(assay(imputed_obj))
#'
#' # Impute using 90% of the row minimum
#' imputed_scaled <- impute_min(se, alpha = 0.9)
#' cat("\nImputed with alpha = 0.9:\n")
#' print(assay(imputed_scaled))
#'
setGeneric("impute_min", function(object, alpha=1) standardGeneric("impute_min"))

#' Impute from a Down-Shifted Normal Distribution
#'
#' @description
#' Performs row-wise imputation by drawing random values from a normal
#' distribution that is shifted to the left and narrower than the distribution of
#' observed values.
#'
#' This method assumes that missing values are primarily from proteins with low
#' abundance (i.e., below the detection limit). The default `shift` and `scale`
#' values are based on those used in the Perseus analysis platform.
#'
#' @details
#' This imputation method should be applied to log-transformed data, as the
#' underlying assumption of a normal distribution is more appropriate in log space.
#' An exception is made for rows with only one or zero observed values, where missing
#' values are imputed with 0.
#'
#' @param object A \code{SummarizedExperiment} object containing data with missing values.
#' @param shift A numeric value specifying how many standard deviations to shift
#'   the mean of the distribution for imputed values. Default is 1.8.
#' @param scale A numeric value to scale the standard deviation of the
#'   distribution for imputed values. Default is 0.3.
#'
#' @return A \code{SummarizedExperiment} object with missing values imputed from a simulated
#'   low-abundance distribution.
#'
#' @export
#' @examples
#' # Create data with NAs, typically representing log-transformed values
#' raw_data <- data.frame(
#'   Gene = c("GENEA", "GENEB"),
#'   SampleA = c(25.1, 28.5),
#'   SampleB = c(25.5, 28.9),
#'   SampleC = c(NA, NA)
#' )
#'
#' se <- create_se(raw_data)
#'
#' # For reproducibility of the random imputation
#' set.seed(123)
#'
#' imputed_obj <- impute_left_dist(se)
#' cat("Data after imputation:\n")
#' print(assay(imputed_obj))
#'
setGeneric("impute_left_dist", function(object, shift = 1.8, scale = 0.3) standardGeneric("impute_left_dist"))

#' Correct for Batch Effects
#'
#' This function uses the `removeBatchEffect` method from the limma package to
#' adjust for known batch effects in the data.
#'
#' @param object A \code{SummarizedExperiment} object. The first assay will be used.
#' @param batch_variable A single character string specifying the column name
#'   in `colData(object)` that contains the batch information.
#' @param bio_variables An optional character vector of column names in
#'   `colData(object)` that correspond to biological variables of interest.
#'   These effects will be preserved during the correction. If `NULL` (the default),
#'   only an intercept is protected.
#'
#' @return A \code{SummarizedExperiment} object with a new assay named
#'   `"corrected"` containing the batch-corrected data.
#'
#' @export
#' @examples
#' # --- Create a sample SummarizedExperiment object ---
#' set.seed(123)
#' counts <- matrix(rnorm(100 * 6, mean = 10, sd = 2), nrow = 100, ncol = 6)
#' colnames(counts) <- paste0("Sample", 1:6)
#'
#' sample_info <- data.frame(
#'   row.names = colnames(counts),
#'   Condition = rep(c("A", "B"), each = 3),
#'   Batch = rep(c("Batch1", "Batch2"), times = 3)
#' )
#'
#' se <- SummarizedExperiment(assays = list(log_intensities = counts),
#'                            colData = sample_info)
#'
#' # --- Run batch correction, preserving the "Condition" variable ---
#' corrected_se <- batch_correct(se,
#'                               batch_variable = "Batch",
#'                               bio_variables = "Condition")
#'
#' # --- Inspect the result ---
#' assayNames(corrected_se) # Now includes "corrected"
#'
setGeneric("batch_correct", function(object, batch_variable, bio_variables = NULL) {
  standardGeneric("batch_correct")
})

## Clustering ############################################################

#' Calculate Principal Components Analysis (PCA)
#'
#' @description
#' Performs a principal component analysis (PCA) on the samples within a
#' `SummarizedExperiment` object.
#'
#' @details
#' This function expects data that has already been cleaned and imputed. Missing
#' values (`NA`) will cause the PCA to fail. It is also highly
#' recommended to perform log-transformation and normalization before PCA.
#' The function will automatically remove proteins (rows) with zero variance
#' across samples before running the analysis.
#'
#' @param object A `SummarizedExperiment` object. The first assay is used.
#' @param condition A character string specifying the column name in
#'   `colData(object)` to use for grouping the samples in the output. If `NA`
#'   (the default), groups are inferred from sample names by removing numeric suffixes.
#'
#' @return
#' A list with two elements:
#' \describe{
#'   \item{`summary`}{A data.table with the standard deviation, proportion of
#'     variance, and cumulative proportion for each component.}
#'   \item{`components`}{A data frame with the first 5 PCA scores (PC1-PC5)
#'     for each sample, along with their assigned condition.}
#' }
#'
#' @export
#' @examples
#' # --- Create a sample SummarizedExperiment object ---
#' set.seed(123)
#' counts <- matrix(rnorm(100 * 6, mean = 10, sd = 2), nrow = 100, ncol = 6)
#' colnames(counts) <- paste0("Sample_", rep(c("A","B"), each=3), "_", 1:3)
#' sample_info <- data.frame(
#'   row.names = colnames(counts),
#'   Group = rep(c("A", "B"), each = 3)
#' )
#' se <- SummarizedExperiment(assays = list(log_intensities = counts),
#'                            colData = sample_info)
#'
#' # --- Run PCA, specifying the group variable ---
#' pca_results <- get_PCs(se, condition = "Group")
#' head(pca_results$components)
#'
setGeneric("get_PCs",
           function(object, condition = NA) {
             standardGeneric("get_PCs")
           }
)

#' Plot Principal Component Analysis Results
#'
#' @description
#' Generates a scatter plot of principal components to visualize the relationships
#' between samples. This function serves as a convenient wrapper that first calls
#' `calculate_pca` and then plots the results using `ggplot2`.
#'
#' @details
#' As this function calls `calculate_pca` internally, the data in the object
#' must be imputed first. It is also recommended to use log-transformed and
#' normalized data for the best results.
#'
#' @param object A `SummarizedExperiment` object.
#' @param condition A character string specifying the column name in the
#'   `condition` slot to use for coloring the points. If `NA` (the default),
#'   it will attempt to guess groups from sample names.
#' @param pc_x A character string for the principal component to plot on the
#'   x-axis (e.g., `"PC1"`). Defaults to `"PC1"`.
#' @param pc_y A character string for the principal component to plot on the
#'   y-axis (e.g., `"PC2"`). Defaults to `"PC2"`.
#'
#' @return A `ggplot` object, which can be further customized or printed.
#'
#' @export
#'
#' @examples
#' # Create a sample SummarizedExperiment object with missing data
#' raw_data <- data.frame(
#'   Gene = c("GENEA", "GENEB", "GENEC", "GENED", "GENEE", "GENEF"),
#'   Control_1 = c(10, 11, 12, 13, 14, 15),
#'   Control_2 = c(10.5, 11.5, 12.5, NA, 14.5, 15.5),
#'   Treatment_1 = c(15, 16, 17, 18, 19, 20),
#'   Treatment_2 = c(15.5, 16.5, 17.5, 18.5, 19.5, 20.5)
#' )
#' cond_df <- data.frame(
#'    SampleID = c("Control_1", "Control_2", "Treatment_1", "Treatment_2"),
#'    group = c("Control", "Control", "Treatment", "Treatment")
#' )
#' se <- create_se(raw_data, sample_metadata = cond_df)
#'
#' # Impute missing values before plotting
#' se_imputed <- impute(se, value = 13.5)
#'
#' # Generate the plot of PC1 vs PC2
#' p1 <- plot_PCs(se_imputed, condition = "group")
#' if (interactive()) {
#'   print(p1)
#' }
#'
#' # Generate a plot of PC1 vs PC3
#' p2 <- plot_PCs(se_imputed, condition = "group", pc_x = "PC1", pc_y = "PC3")
#' if (interactive()) {
#'   print(p2)
#' }
#'
setGeneric("plot_PCs",
           function(object, condition = NA, pc_x = "PC1", pc_y = "PC2") {
             standardGeneric("plot_PCs")
           }
)

#' Plot a Hierarchical Clustering Dendrogram of Samples
#'
#' @description
#' Performs hierarchical clustering on the samples based on their protein
#' abundance profiles and generates a dendrogram plot using the `ggdendro`
#' package.
#'
#' @details
#' This function expects clean, imputed data. Missing values (`NA`) will cause
#' an error. For meaningful biological results, it is highly recommended to use
#' data that has been log-transformed and normalized before clustering.
#'
#' @param object A `SummarizedExperiment` object. The data should be imputed.
#' @param dist_method The distance measure to be used by `stats::dist`. Common
#'   options include `"euclidean"`, `"maximum"`, `"manhattan"`. Defaults to `"euclidean"`.
#' @param hclust_method The agglomeration method to be used by `stats::hclust`.
#'   Common options include `"complete"`, `"ward.D2"`, `"average"`. Defaults to `"complete"`.
#'
#' @return A `ggplot` object representing the dendrogram, which can be further
#'   customized.
#'
#' @export
#'
#' @examples
#' # Create a sample SummarizedExperiment object
#' raw_data <- data.frame(
#'   Gene = c("GENEA", "GENEB", "GENEC", "GENED"),
#'   SampleA = c(10, 20, 15, 12),
#'   SampleB = c(11, 21, 16, 13), # Similar to A
#'   SampleC = c(25, 10, 30, 5),
#'   SampleD = c(26, 11, 31, 6)  # Similar to C
#' )
#' se <- create_se(raw_data)
#'
#' # Run with default methods. We expect A/B and C/D to cluster together.
#' p1 <- plot_hierarchical_cluster(se)
#' if (interactive()) {
#'   print(p1)
#' }
#'
#' # Run with different methods
#' p2 <- plot_hierarchical_cluster(se, dist_method = "manhattan", hclust_method = "ward.D2")
#' if (interactive()) {
#'   print(p2)
#' }
#'
setGeneric("plot_hierarchical_cluster",
           function(object, dist_method = "euclidean", hclust_method = "complete") {
             standardGeneric("plot_hierarchical_cluster")
           }
)

#' Calculate UMAP Dimensionality Reduction
#'
#' @description
#' Performs Uniform Manifold Approximation and Projection (UMAP) on the samples
#' to generate a 2D embedding. This is a powerful non-linear method for
#' visualizing sample relationships. This function serves as a wrapper for the
#' `umap::umap` function.
#'
#' @details
#' This function expects clean, imputed data. Missing values (`NA`) will cause
#' an error. For meaningful results, it is highly recommended to use data that has
#' been log-transformed and normalized.
#'
#' **Important:** UMAP is a stochastic algorithm, meaning it will produce slightly
#' different results each time it is run. For reproducible results, you **must set
#' a seed** (e.g., `set.seed(123)`) before calling this function.
#'
#' This function requires the `umap` package to be installed from CRAN.
#'
#' @param object A `SummarizedExperiment` object. The data should be imputed.
#' @param condition A character string specifying the column name in the
#'   `condition` slot to use for labeling points. If `NULL` (the default), it
#'   will attempt to guess groups from sample names.
#' @param neighbors The size of the local neighborhood UMAP will look at. This is a
#'   key hyperparameter affecting the balance between local and global structure.
#'   Defaults to 15.
#' @param ... Additional arguments passed on to the `umap::umap` function (e.g.,
#'   `min_dist`, `n_epochs`, `metric`).
#'
#' @return A `data.table` with columns for `Sample`, `UMAP1`, `UMAP2`, and `Condition`.
#'
#' @export
#'
setGeneric("get_umap",
           function(object, condition = NA, neighbors = 15, ...) {
             standardGeneric("get_umap")
           }
)

#' Plot UMAP Dimensionality Reduction Results
#'
#' @description
#' Generates a 2D scatter plot of UMAP results to visualize sample relationships.
#' This is a convenient wrapper function that first calls `get_umap` to calculate
#' the coordinates and then plots them using `ggplot2`.
#'
#' @details
#' This function requires imputed, log-transformed, and normalized data for the
#' best results, as these are the best inputs for the wrapped `get_umap` function.
#'
#' **Important:** UMAP is a stochastic algorithm. For a reproducible plot, you
#' **must set a seed** (e.g., `set.seed(123)`) *before* calling this function.
#'
#' This function requires the `umap` and `ggplot2` packages.
#'
#' @param object A `SummarizedExperiment` object. The data should be imputed.
#' @param condition A character string specifying the column name in the
#'   `condition` slot to use for coloring the points. If `NULL` (the default),
#'   `get_umap` will attempt to guess groups from sample names.
#' @param neighbors The size of the local neighborhood for the UMAP calculation.
#' @param ... Additional arguments passed on to `get_umap`, and in turn to
#'   the `umap::umap` function (e.g., `min_dist`, `metric`).
#'
#' @return A `ggplot` object, which can be further customized or printed.
#'
#' @export
#'
#' @examples
#' # This example requires the 'umap' package
#' if (requireNamespace("umap", quietly = TRUE)) {
#'   # Create a sample SummarizedExperiment object
#'   raw_data <- data.frame(
#'     Gene = paste0("GENE", 1:10),
#'     Control_1 = rnorm(10, 10), Control_2 = rnorm(10, 10),
#'     Treat_A_1 = rnorm(10, 12), Treat_A_2 = rnorm(10, 12),
#'     Treat_B_1 = rnorm(10, 15), Treat_B_2 = rnorm(10, 15)
#'   )
#'   se <- create_se(raw_data)
#'
#'   # For reproducible UMAP results, set a seed!
#'   set.seed(42)
#'
#'   # Generate the UMAP plot
#'   p <- plot_umap(se, neighbors = 3)
#'
#'   # The plot can be printed in an interactive session
#'   if (interactive()) {
#'     print(p)
#'   }
#' }
#'
setGeneric("plot_umap",
           function(object, condition = NA, neighbors = 15, ...) {
             standardGeneric("plot_umap")
           }
)

## Differential Expression ###############################################

## Heatmap ############################################################

#' Plot a Proteomics Heatmap
#'
#' @description
#' Generates a heatmap of protein expression. The function can operate in two modes
#' depending on the `condition` argument.
#'
#' 1.  **Individual Samples (default):** If `condition` is `NULL`, a heatmap of
#'     all individual samples is generated, with samples clustered by similarity.
#' 2.  **Summarized Conditions:** If a `condition` is provided, the function first
#'     calculates the mean expression for each protein across the replicates
#'     within each condition group, then generates a heatmap of these mean values.
#'
#' In both modes, the data is row-wise Z-scored before plotting.
#'
#' @param object A `SummarizedExperiment` object. It is recommended to use data that has been
#'   log-transformed and imputed.
#' @param protmeta_col A character string specifying the name of the column in
#'   the `rowData` slot to use for protein labels on the heatmap rows.
#'   Defaults to the first column of `rowData`.
#' @param genes Optional. A character vector of gene or protein names to subset
#'   the data to. Only these proteins will be included in the heatmap. The match
#'   is case-insensitive.
#' @param title Optional. A character string for the plot title.
#' @param condition Optional. A character string specifying a column name in the
#'   `colData` slot. If provided, triggers the summarized heatmap mode.
#' @param cluster_cols Logical indicating whether heatmap columns should be clustered.
#' @param cluster_rows Logical indicating whether heatmap rows should be clustered.
#' @return A `ggplot` object representing the heatmap.
#'
#' @export
#'
setGeneric("plot_proteomics_heatmap",
           function(object, protmeta_col = NULL, genes = NULL, title = NULL, condition = NULL,
                    cluster_cols = FALSE, cluster_rows = FALSE) {
             standardGeneric("plot_proteomics_heatmap")
           }
)

## Protein Comparison ############################################################
