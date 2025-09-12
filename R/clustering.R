#' Calculate Principal Components Analysis (PCA)
#'
#' @description
#' Performs a principal component analysis (PCA) on the samples within a
#' `ProtData` object. This is a common method for quality control and for
#' visualizing the relationships between samples.
#'
#' @details
#' This function expects data that has already been cleaned and imputed. Missing
#' values (`NA`) in the data will cause the PCA to fail. It is also highly
#' recommended to perform log-transformation and normalization before PCA.
#'
#' The function will automatically remove proteins (rows) with no variation
#' across samples before running the analysis.
#'
#' @param object A `ProtData` object. The data should be imputed and ideally
#'   normalized and log-transformed.
#' @param condition A character string specifying the column name in the
#'   `condition` slot to use for grouping the samples in the output. If `NA`
#'   (the default), it will attempt to guess groups by removing numeric suffixes
#'   from sample names.
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
#'
setGeneric("get_PCs",
           function(object, condition = NA) {
             standardGeneric("get_PCs")
           }
)

setMethod("get_PCs", "ProtData",
          function(object, condition) {

            # --- 1. Data Preparation (from refactored code) ---
            data_for_pca <- object@data
            data_for_pca <- data_for_pca[rowSums(abs(data_for_pca), na.rm = TRUE) > 0, ]

            if (any(is.na(data_for_pca))) {
              stop("Missing values (NA) found in data. Please impute before running PCA.")
            }

            pca_data_filtered <- t(data_for_pca)
            variances <- apply(pca_data_filtered, 2, var, na.rm = TRUE)
            pca_data_filtered <- pca_data_filtered[, variances > 1e-9]

            if (ncol(pca_data_filtered) < 2) {
              stop("Not enough features with variance to perform PCA.")
            }

            # --- 2. PCA and Output Formatting (Your Original Logic) ---
            out <- list()
            pca <- stats::prcomp(pca_data_filtered, center = TRUE, scale. = TRUE)
            out$summary <- data.table::as.data.table(t(summary(pca)$importance), keep.rownames = TRUE)
            data.table::setnames(out$summary, c('component','stdv','percent','cumulative'))
            out$summary$percent <- round(out$summary$percent * 100, digits = 2)

            pca_df <- as.data.frame(pca$x)[, 1:5]
            pca_df$Sample <- rownames(pca_df)

            if (is.na(condition)) {
              pca_df$Condition <- gsub('_[0-9]+$', '', rownames(pca_df))
            } else {
              if (!condition %in% names(object@condition)) {
                stop("'", condition, "' not found in the condition slot.")
              }
              pca_df$Condition <- object@condition[[condition]]
            }

            out$components <- pca_df
            return(out)
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
#' @param object A `ProtData` object.
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
#' # Create a sample ProtData object with missing data
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
#' pd_obj <- create_protdata(dat = raw_data, condition = cond_df)
#'
#' # Impute missing values before plotting
#' pd_obj_imputed <- impute(pd_obj, value = 13.5)
#'
#' # Generate the plot of PC1 vs PC2
#' p1 <- plot_pca(pd_obj_imputed, condition = "group")
#' if (interactive()) {
#'   print(p1)
#' }
#'
#' # Generate a plot of PC1 vs PC3
#' p2 <- plot_pca(pd_obj_imputed, condition = "group", pc_x = "PC1", pc_y = "PC3")
#' if (interactive()) {
#'   print(p2)
#' }
#'
setGeneric("plot_PCs",
           function(object, condition = NA, pc_x = "PC1", pc_y = "PC2") {
             standardGeneric("plot_PCs")
           }
)

setMethod("plot_PCs", "ProtData",
          function(object, condition, pc_x, pc_y) {

            # Step 1: Call the calculation method to get PCA results
            pca_results <- get_PCs(object, condition = condition)

            pca_df <- pca_results$components
            pca_summary <- pca_results$summary

            # Check if the requested PCs are available
            if (!pc_x %in% names(pca_df)) stop("X-axis component not found:", pc_x)
            if (!pc_y %in% names(pca_df)) stop("Y-axis component not found:", pc_y)

            # Step 2: Dynamically create axis labels with variance explained
            percent_x <- pca_summary[pca_summary$component == pc_x, ]$percent
            percent_y <- pca_summary[pca_summary$component == pc_y, ]$percent

            xlab_text <- paste0(pc_x, " (", percent_x, "%)")
            ylab_text <- paste0(pc_y, " (", percent_y, "%)")

            # Step 3: Create the plot
            # Use the .data pronoun for tidy evaluation of string variable names
            p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = .data[[pc_x]], y = .data[[pc_y]], color = Condition)) +
              ggplot2::geom_point(size = 4, alpha = 0.8) +
              ggplot2::xlab(xlab_text) +
              ggplot2::ylab(ylab_text) +
              ggplot2::theme_classic() +
              ggplot2::labs(
                title = "PCA of Samples",
                color = if(!is.na(condition)) condition else "Condition"
              )

            return(p)
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
#' @param object A `ProtData` object. The data should be imputed.
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
#' # Create a sample ProtData object
#' raw_data <- data.frame(
#'   Gene = c("GENEA", "GENEB", "GENEC", "GENED"),
#'   SampleA = c(10, 20, 15, 12),
#'   SampleB = c(11, 21, 16, 13), # Similar to A
#'   SampleC = c(25, 10, 30, 5),
#'   SampleD = c(26, 11, 31, 6)  # Similar to C
#' )
#' pd_obj <- create_protdata(dat = raw_data)
#'
#' # Run with default methods. We expect A/B and C/D to cluster together.
#' p1 <- plot_hierarchical_cluster(pd_obj)
#' if (interactive()) {
#'   print(p1)
#' }
#'
#' # Run with different methods
#' p2 <- plot_hierarchical_cluster(pd_obj, dist_method = "manhattan", hclust_method = "ward.D2")
#' if (interactive()) {
#'   print(p2)
#' }
#'
setGeneric("plot_hierarchical_cluster",
           function(object, dist_method = "euclidean", hclust_method = "complete") {
             standardGeneric("plot_hierarchical_cluster")
           }
)

setMethod("plot_hierarchical_cluster", "ProtData",
          function(object, dist_method, hclust_method) {

            # --- 1. Data Preparation ---
            cluster_data <- object@data

            # Check for missing values; this is critical.
            if (any(is.na(cluster_data))) {
              stop("Missing values (NA) found. Please impute data before clustering.")
            }

            # Filter rows that are all zero to avoid issues with distance calculations
            cluster_data <- cluster_data[rowSums(abs(cluster_data)) > 0, ]

            # --- 2. Perform Clustering ---

            # Calculate distance matrix on samples (columns need to be transposed)
            dist_mat <- stats::dist(t(cluster_data), method = dist_method)

            # Perform hierarchical clustering
            hc_cluster <- stats::hclust(dist_mat, method = hclust_method)

            # --- 3. Generate Plot ---
            p <- ggdendro::ggdendrogram(hc_cluster, rotate = TRUE) +
              ggplot2::labs(
                title = "Hierarchical Clustering of Samples",
                y = "Distance",
                x = "Sample"
              )

            return(p)
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
#' @param object A `ProtData` object. The data should be imputed.
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

setMethod("get_umap", "ProtData",
          function(object, condition, neighbors, ...) {
            ttt<<- condition
            # --- 1. Input Validation and Data Prep ---
            if (!requireNamespace("umap", quietly = TRUE)) stop("Please install the 'umap' package.")

            umap_data <- object@data

            if (any(is.na(umap_data))) {
              stop("Missing values (NA) found. Please impute data before running UMAP.")
            }

            umap_data <- umap_data[rowSums(abs(umap_data), na.rm = TRUE) > 0, ]

            # --- 2. Run UMAP ---
            # The function does NOT set a seed; the user must do this for reproducibility.
            umap_result <- umap::umap(t(umap_data), n_neighbors = neighbors, ...)

            # --- 3. Format Output ---
            umap_df <- data.table::as.data.table(umap_result$layout, keep.rownames = TRUE)
            data.table::setnames(umap_df, c('Sample', 'UMAP1', 'UMAP2'))

            if (is.na(condition)) {
              umap_df$Condition <- gsub('_[0-9]+$', '', umap_df$Sample)
            } else {
              if (!condition %in% names(object@condition)) {
                stop("'", condition, "' not found in the condition slot.")
              }
              # Safely get condition data, ensuring row order matches
              umap_df$Condition <- object@condition[[condition]]
            }

            return(umap_df)
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
#' @param object A `ProtData` object. The data should be imputed.
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
#'   # Create a sample ProtData object
#'   raw_data <- data.frame(
#'     Gene = paste0("GENE", 1:10),
#'     Control_1 = rnorm(10, 10), Control_2 = rnorm(10, 10),
#'     Treat_A_1 = rnorm(10, 12), Treat_A_2 = rnorm(10, 12),
#'     Treat_B_1 = rnorm(10, 15), Treat_B_2 = rnorm(10, 15)
#'   )
#'   pd_obj <- create_protdata(dat = raw_data)
#'
#'   # For reproducible UMAP results, set a seed!
#'   set.seed(42)
#'
#'   # Generate the UMAP plot
#'   p <- plot_umap(pd_obj, n_neighbors = 3)
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

setMethod("plot_umap", "ProtData",
          function(object, condition, neighbors, ...) {

            # Step 1: Call the calculation method to get UMAP coordinates.
            # Pass the ellipsis (...) to allow for more umap arguments.
            umap_df <- get_umap(
              object = object,
              condition = condition,
              n_neighbors = neighbors,
              ...
            )

            # Step 2: Create the plot
            p <- ggplot2::ggplot(umap_df, ggplot2::aes(x = UMAP1, y = UMAP2, color = Condition)) +
              ggplot2::geom_point(size = 4, alpha = 0.8) +
              ggplot2::theme_classic() +
              ggplot2::labs(
                title = "UMAP of Samples",
                color = if(!is.na(condition)) condition else "Condition"
              )

            return(p)
          }
)










#plsda



# --- 1. Installation of the mixOmics package ---
# mixOmics is a Bioconductor package.
if (!requireNamespace("mixOmics", quietly = TRUE)) {
  if(!require("BiocManager")){
    install.packages("BiocManager")
  }
  BiocManager::install("mixOmics")
}


# --- 2. Generic Function Definition ---
#' Generate a PLS-DA scores plot for a ProtData object
#'
#' @param object The ProtData object.
#' @param group_variable A string naming the column in the condition slot that contains the group labels.
#' @param n_components The number of components for the PLS-DA model.
setGeneric("plot_plsda",
           def = function(object, group_variable, n_components = 2) {
             standardGeneric("plot_plsda")
           }
)


# --- 3. S4 Method Implementation ---
#' @describeIn plot_plsda PLS-DA method for ProtData class.
setMethod("plot_plsda",
          signature(object = "ProtData", group_variable = "character"),
          function(object, group_variable, n_components = 2) {

            # --- Input Validation ---
            if (!group_variable %in% names(object@condition)) {
              stop("'", group_variable, "' not found in the condition slot of the object.")
            }
            if (any(is.na(object@data))) {
              stop("Missing values (NA) detected in the data slot. Please impute or filter data before running PLS-DA.")
            }

            groups <- as.factor(object@condition[[group_variable]])
            if (length(levels(groups)) < 2) {
              stop("PLS-DA requires at least two groups in your '", group_variable, "' column.")
            }

            cat("Running PLS-DA for group variable:", group_variable, "\n")

            # --- Prepare Data for mixOmics ---
            # mixOmics requires samples as rows and features (proteins) as columns.
            # We must transpose the data matrix.
            X <- t(object@data)
            Y <- groups

            # --- Run PLS-DA ---
            plsda_result <- mixOmics::plsda(X, Y, ncomp = n_components)

            # --- Generate the Plot ---
            # plotIndiv is mixOmics's powerful plotting function.
            # We use style='ggplot2' to get a customizable ggplot object back.
            scores_plot <- mixOmics::plotIndiv(
              plsda_result,
              group = Y,
              ind.names = FALSE, # We don't need sample names on the plot
              legend = TRUE,
              ellipse = TRUE,    # Draw confidence ellipses around the groups
              style = 'ggplot2',
              title = paste("PLS-DA Scores Plot -", Sys.Date())
            )

            return(scores_plot)
          }
)
