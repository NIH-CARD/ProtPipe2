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
plot_proteomics_heatmap <- function(PD, protmeta_col = NULL, genes = NULL, title = NULL, condition = NULL) {
  # Check input class
  if (!inherits(PD, "ProtData")) {
    stop("Error: 'PD' must be of class ProtData")
  }

  #Extract data
  intensities <- PD@data
  prot_meta <- PD@prot_meta
  condition_file <- PD@condition

  if(is.null(protmeta_col)){
    protmeta_col <- names(prot_meta)[1]
  }

  if(is.null(title)){
    title <- "Proteomics Heatmap"
  }

  # Check if protmeta_col exists
  if (!protmeta_col %in% colnames(prot_meta)) {
    stop(paste("Column", protmeta_col, "not found in PD@prot_meta"))
  }

  # Optionally subset by genes
  if (!is.null(genes)) {
    match_genes <- tolower(prot_meta[[protmeta_col]]) %in% tolower(genes)
    intensities <- intensities[match_genes, , drop = FALSE]
    prot_meta <- prot_meta[match_genes, , drop = FALSE]
  }

  # Impute NA/NaN/Inf with 0
  intensities <- as.matrix(intensities)
  #intensities[!is.finite(intensities)] <- 0

  # Make unique row labels using protmeta_col
  labels <- make.unique(as.character(prot_meta[[protmeta_col]]))
  rownames(intensities) <- labels

  show_rownames <- nrow(intensities) < 100

  if(is.null(condition)){
    col_labels <- names(intensities)
  }else{
    if (condition %in% colnames(condition_file)){
      col_labels <- condition_file[[condition]]
    }else{
      stop("the selected condition does not appear in the condition file")
    }
  }

  pheatmap::pheatmap(
    mat = intensities,
    scale = 'row',
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = show_rownames,
    show_colnames = TRUE,
    labels_col = col_labels,
    fontsize_row = 8,
    fontsize_col = 10,
    main = title
  )
}
