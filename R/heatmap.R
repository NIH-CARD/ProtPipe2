#' Title
#'
#' @param PD
#' @param protmeta_col
#' @param genes
#'
#' @return
#' @export
#'
#' @examples
# plot_proteomics_heatmap <- function(PD, protmeta_col = NULL, genes = NULL, title = NULL, condition = NULL) {
#   # Check input class
#   if (!inherits(PD, "ProtData")) {
#     stop("Error: 'PD' must be of class ProtData")
#   }
#
#   #Extract data
#   intensities <- PD@data
#   prot_meta <- PD@prot_meta
#   condition_file <- PD@condition
#
#   if(is.null(protmeta_col)){
#     protmeta_col <- names(prot_meta)[1]
#   }
#
#   if(is.null(title)){
#     title <- "Proteomics Heatmap"
#   }
#
#   # Check if protmeta_col exists
#   if (!protmeta_col %in% colnames(prot_meta)) {
#     stop(paste("Column", protmeta_col, "not found in PD@prot_meta"))
#   }
#
#   # Optionally subset by genes
#   if (!is.null(genes)) {
#     match_genes <- tolower(prot_meta[[protmeta_col]]) %in% tolower(genes)
#     intensities <- intensities[match_genes, , drop = FALSE]
#     prot_meta <- prot_meta[match_genes, , drop = FALSE]
#   }
#
#   # Impute NA/NaN/Inf with 0
#   intensities <- as.matrix(intensities)
#   #intensities[!is.finite(intensities)] <- 0
#
#   # Make unique row labels using protmeta_col
#   labels <- make.unique(as.character(prot_meta[[protmeta_col]]))
#   rownames(intensities) <- labels
#
#   show_rownames <- nrow(intensities) < 100
#
#   if(is.null(condition)){
#     col_labels <- names(intensities)
#   }else{
#     if (condition %in% colnames(condition_file)){
#       col_labels <- condition_file[[condition]]
#     }else{
#       stop("the selected condition does not appear in the condition file")
#     }
#   }
#
#   pheatmap::pheatmap(
#     mat = intensities,
#     scale = 'row',
#     cluster_rows = FALSE,
#     cluster_cols = TRUE,
#     show_rownames = show_rownames,
#     show_colnames = TRUE,
#     labels_col = col_labels,
#     fontsize_row = 8,
#     fontsize_col = 10,
#     main = title
#   )
# }
#
# # Make sure you have the necessary packages installed
# # install.packages(c("ggplot2", "tidyr", "dplyr", "tibble"))

plot_proteomics_heatmap <- function(PD, protmeta_col = NULL, genes = NULL, title = NULL, condition = NULL) {
  #PD <- log2_transform(PD)
  # Check input class
  if (!inherits(PD, "ProtData")) {
    stop("Error: 'PD' must be of class ProtData")
  }

  # --- 1. Data Extraction and Filtering (same as before) ---
  intensities <- PD@data
  prot_meta <- PD@prot_meta
  condition_file <- PD@condition

  if (is.null(protmeta_col)) {
    protmeta_col <- names(prot_meta)[1]
  }

  if (is.null(title)) {
    title <- "Proteomics Heatmap"
  }

  if (!protmeta_col %in% colnames(prot_meta)) {
    stop(paste("Column", protmeta_col, "not found in PD@prot_meta"))
  }

  if (!is.null(genes)) {
    match_genes <- tolower(prot_meta[[protmeta_col]]) %in% tolower(genes)
    intensities <- intensities[match_genes, , drop = FALSE]
    prot_meta <- prot_meta[match_genes, , drop = FALSE]
  }

  if (nrow(intensities) == 0) {
    stop("No data remains after filtering. Cannot generate heatmap.")
  }

  # --- 2. Data Preparation for ggplot2 ---

  # Scale data by row (protein) to get Z-scores
  # The t() transposes the matrix, scale() works on columns, then t() transposes back
  scaled_intensities <- t(base::scale(t(as.matrix(intensities))))

  # Replace any NA/NaN/Inf from scaling (e.g., from rows with 0 variance) with 0
  scaled_intensities[!is.finite(scaled_intensities)] <- 0

  # Perform hierarchical clustering on the columns (samples) to get their order
  col_clustering <- hclust(dist(t(scaled_intensities)))
  col_order <- colnames(scaled_intensities)[col_clustering$order]

  # Get protein labels and ensure they are unique
  protein_labels <- make.unique(as.character(prot_meta[[protmeta_col]]))

  # Convert the scaled matrix to a long-format data frame (tibble)
  heatmap_data <- as.data.frame(scaled_intensities) %>%
    tibble::rownames_to_column(var = "Protein") %>%
    dplyr::mutate(Protein = factor(protein_labels, levels = rev(protein_labels))) %>% # Use labels, reverse order for ggplot
    tidyr::pivot_longer(
      cols = -Protein,
      names_to = "Sample",
      values_to = "ZScore"
    ) %>%
    dplyr::mutate(Sample = factor(Sample, levels = col_order)) # Apply clustering order

  heatmap_data<<- heatmap_data
  # --- 3. Handle Column Labels from Condition File ---

  # Default to sample names if no condition is provided
  plot_col_labels <- levels(heatmap_data$Sample)

  if (!is.null(condition)) {
    if (condition %in% colnames(condition_file)) {
      # Create a named vector for mapping sample names to condition labels
      # Assumes condition_file has a column with sample names that match colnames(intensities)
      # Let's assume the first column of condition_file contains the sample names
      sample_name_col <- names(condition_file)[1]
      label_map <- setNames(condition_file[[condition]], condition_file[[sample_name_col]])
      plot_col_labels <- label_map[levels(heatmap_data$Sample)]
    } else {
      stop("The selected condition does not appear in the condition file")
    }
  }

  # --- 4. Create the ggplot Heatmap ---
  max_val <- ceiling(max(heatmap_data$ZScore))
  min_val <- floor(min(heatmap_data$ZScore))
 g <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = Sample, y = Protein, fill = ZScore)) +
  ggplot2::geom_tile(color = "black", size = 0.0) +
  ggplot2::theme_classic() +
  ggplot2::scale_fill_gradient2(
    low = "skyblue", high = "tomato1", mid = "beige",
    midpoint = 0,
    #space = "Lab", breaks = c(-max_abs_val, max_abs_val),
    name = "protein\nz-score",
    limits = c(min_val, max_val),
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
    axis.text.y = ggplot2::element_text(size = 10),
    axis.title = ggplot2::element_blank()
  )

  # Conditionally hide row labels if there are too many
  if (nrow(intensities) >= 100) {
    g <- g + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                            axis.ticks.y = ggplot2::element_blank())
  }

  return(g)
}


