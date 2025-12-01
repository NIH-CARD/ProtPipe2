#' @describeIn get_PCs Method for SummarizedExperiment objects.
#' @export
setMethod("get_PCs", "SummarizedExperiment",
          function(object, condition) {

            # --- 1. Data Preparation ---
            # Use the first assay by default
            data_for_pca <- assay(object)

            # Stop if data contains missing values
            if (anyNA(data_for_pca)) {
              stop("Missing values (NA) found. Please impute before running PCA.")
            }

            # Transpose the data so that samples are rows and proteins are columns
            pca_data_transposed <- t(data_for_pca)

            # Filter out proteins (now columns) with zero variance
            variances <- apply(pca_data_transposed, 2, var)
            pca_data_filtered <- pca_data_transposed[, variances > 1e-9]

            if (ncol(pca_data_filtered) < 2) {
              stop("Not enough features with variance to perform PCA.")
            }

            # --- 2. PCA and Output Formatting ---
            pca <- stats::prcomp(pca_data_filtered, center = TRUE, scale. = TRUE)

            # Format the variance summary table
            pca_summary <- data.table::as.data.table(t(summary(pca)$importance), keep.rownames = TRUE)
            data.table::setnames(pca_summary, c('component','stdv','percent','cumulative'))
            pca_summary$percent <- round(pca_summary$percent * 100, digits = 2)

            # Format the principal components data frame
            pca_df <- as.data.frame(pca$x)[, 1:5]
            pca_df$Sample <- rownames(pca_df)

            # Add the condition information
            if (is.na(condition)) {
              pca_df$Condition <- gsub('_[0-9]+$', '', rownames(pca_df))
            } else {
              if (!condition %in% colnames(colData(object))) {
                stop("'", condition, "' not found in colData(object).")
              }
              # Match condition info to the PCA results
              pca_df$Condition <- colData(object)[pca_df$Sample, condition]
            }

            # --- 3. Return the List ---
            out <- list(
              summary = pca_summary,
              components = pca_df
            )

            # Note: Logging is omitted here as we are not returning the SE object
            return(out)
          }
)



#' @describeIn plot_PCs Method for SummarizedExperiment objects.
#' @export
setMethod("plot_PCs", "SummarizedExperiment",
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

#' @describeIn plot_hierarchical_cluster Method for SummarizedExperiment objects.
#' @export
setMethod("plot_hierarchical_cluster", "SummarizedExperiment",
          function(object, dist_method, hclust_method) {

            # --- 1. Data Preparation ---
            cluster_data <- assay(object)

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

#' @describeIn get_umap Method for SummarizedExperiment objects.
#' @export
setMethod("get_umap", "SummarizedExperiment",
          function(object, condition, neighbors) {
            # --- 1. Input Validation and Data Prep ---
            if (!requireNamespace("umap", quietly = TRUE)) stop("Please install the 'umap' package.")

            umap_data <- assay(object)

            if (any(is.na(umap_data))) {
              stop("Missing values (NA) found. Please impute data before running UMAP.")
            }

            umap_data <- umap_data[rowSums(abs(umap_data), na.rm = TRUE) > 0, ]

            # --- 2. Run UMAP ---
            # The function does NOT set a seed; the user must do this for reproducibility.
            umap_result <- umap::umap(t(umap_data), n_neighbors = neighbors)

            # --- 3. Format Output ---
            umap_df <- data.table::as.data.table(umap_result$layout, keep.rownames = TRUE)
            data.table::setnames(umap_df, c('Sample', 'UMAP1', 'UMAP2'))

            if (is.na(condition)) {
              umap_df$Condition <- gsub('_[0-9]+$', '', umap_df$Sample)
            } else {
              if (!condition %in% names(colData(object))) {
                stop("'", condition, "' not found in the condition slot.")
              }
              # Safely get condition data, ensuring row order matches
              umap_df$Condition <- colData(object)[[condition]]
            }

            return(umap_df)
          }
)

#' @describeIn plot_umap Method for SummarizedExperiment objects.
#' @export
setMethod("plot_umap", "SummarizedExperiment",
          function(object, condition, neighbors) {

            # Step 1: Call the calculation method to get UMAP coordinates.
            # Pass the ellipsis (...) to allow for more umap arguments.
            umap_df <- get_umap(
              object = object,
              condition = condition,
              neighbors = neighbors
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










#plsda #########################################



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
