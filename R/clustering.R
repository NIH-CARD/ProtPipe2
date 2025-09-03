#' @importFrom magrittr %>%

## PCA
#' Title
#'
#' @param PD
#'
#' @return
#' @export
#'
#' @examples
get_PCs <- function(PD, condition = NA) {
  PD <- ProtPipe::impute(PD, 0) #%>% ProtPipe::log2_transform()
  log2_cluster_data <- getData(PD) %>%
    dplyr::filter(rowSums(abs(.), na.rm = TRUE) > 0)
  log2_cluster_data <<- log2_cluster_data
  out <- list()

  ##PCA and plot
  pca_data <- t(log2_cluster_data)

  # remove zero variance columns
  variances <- apply(pca_data, 2, var)
  pca_data_filtered <- pca_data[, variances > 0]

  pca_data <<- pca_data
  pca=stats::prcomp(pca_data_filtered, center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
  out$summary <- data.table::as.data.table(t(summary(pca)$importance), keep.rownames=T)
  data.table::setnames(out$summary, c('component','stdv','percent','cumulative'))
  out$summary$percent=round(out$summary$percent*100, digits = 2)
  pca_df = as.data.frame(pca$x)[,1:5]
  pca_df$Sample=rownames(pca_df)
  if(is.na(condition)){
    pca_df$Condition=gsub('_[0-9]+$','',rownames(pca_df))
  }else {
    pca_df$Condition=getCondition(PD)[[condition]]
  }
  out$components <- pca_df
  return(out)
}


#' Title
#'
#' @param PD
#'
#' @return
#' @export
#'
#' @examples
plot_PCs <- function(PD, condition = NA) {
  PCA <- get_PCs(PD, condition)
  p <- ggplot2::ggplot(PCA$components, ggplot2::aes(x = PC1, y = PC2, color = Condition)) +
    ggplot2::geom_point(size=4) +
    ggplot2::xlab(paste0("PC1","(",PCA$summary$percent[1],"%)")) +
    ggplot2::ylab(paste0("PC2","(",PCA$summary$percent[2],"%)")) +
    ggplot2::theme_classic()
  return(p)
}

## HC cluster
#' Title
#'
#' @param PD
#'
#' @return
#' @export
#'
#' @examples
plot_hierarchical_cluster <- function(PD) {
  DT <- getData(PD)
  cluster_data <- DT %>%
    dplyr::select_if(is.numeric) %>%
    # Replace NA values with 0
    dplyr::mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))  %>%
    # Filter rows where the sum of numeric values is greater than 0
    dplyr::filter(rowSums(.) > 0)

  # Apply log2 transformation
  log2_cluster_data <- cluster_data #%>%
    #dplyr::mutate(dplyr::across(dplyr::everything(), ~ log2(. + 1)))

  dist_mat <- stats::dist(t(log2_cluster_data)) #
  hc_cluster <- stats::hclust(dist_mat,method = "complete")
  g <- ggdendro::ggdendrogram(hc_cluster, rotate=TRUE) + ggplot2::labs(title='Hierarchical clustering')
  return(g)
}

## umap
#' Title
#'
#' @param PD
#' @param neighbors
#'
#' @return
#' @export
#'
#' @examples
get_umap <- function(PD, neighbors = 15, condition = NA) {
  DT <- getData(PD)
  ##cluster data(na=0)
  cluster_data <- DT %>%
    dplyr::select_if(is.numeric)%>%
    # Replace NA values with 0
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ ifelse(is.na(.), 0, .)))  %>%
    # Filter rows where the sum of numeric values is greater than 0
    dplyr::filter(rowSums(.) > 0)

  # Apply log2 transformation
  log2_cluster_data <- cluster_data #%>%
    #dplyr::mutate(across(everything(), ~ log2(. + 1)))

  set.seed(100)
  DT.umap <- umap::umap(t(log2_cluster_data), n_neighbors=neighbors)
  DT.out <- data.table::as.data.table(DT.umap$layout, keep.rownames=TRUE)
  data.table::setnames(DT.out, c('Sample', 'UMAP1', 'UMAP2'))
  if(is.na(condition)){
    DT.out$Condition=gsub('_[0-9]+$','',DT.out$Sample)
  }else {
    DT.out$Condition=getCondition(PD)[[condition]]
  }
  return(DT.out[])
}

#' Title
#'
#' @param PD
#'
#' @return
#' @export
#'
#' @examples
plot_umap <- function(PD, neighbors = 15, condition = NA) {
  DT <- get_umap(PD, neighbors = neighbors, condition = condition)
  g <- ggplot2::ggplot(DT, ggplot2::aes(x=UMAP1, y=UMAP2, color=Condition)) +
    ggplot2::geom_point(size=4) +
    ggplot2::theme_classic()
  return(g)
}

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
#' @export
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
