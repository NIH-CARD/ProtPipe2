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
#' @param PD A `ProtData` object. It is recommended to use data that has been
#'   log-transformed and imputed.
#' @param protmeta_col A character string specifying the name of the column in
#'   the `@prot_meta` slot to use for protein labels on the heatmap rows.
#'   Defaults to the first column of `@prot_meta`.
#' @param genes Optional. A character vector of gene or protein names to subset
#'   the data to. Only these proteins will be included in the heatmap. The match
#'   is case-insensitive.
#' @param title Optional. A character string for the plot title.
#' @param condition Optional. A character string specifying a column name in the
#'   `@condition` slot. If provided, triggers the summarized heatmap mode.
#'
#' @return A `ggplot` object representing the heatmap.
#'
#' @export
#'
#' @examples
#' # Create sample data for the constructor
#' raw_data <- data.frame(
#'   ProteinID = paste0("P", 1:3), Gene = c("GENEA", "GENEB", "GENEC"),
#'   Ctrl_1 = c(10, 20, 15), Ctrl_2 = c(12, 22, 14),
#'   Treat_1 = c(25, 10, 5), Treat_2 = c(27, 12, 6)
#' )
#' cond_df <- data.frame(
#'    SampleID = c("Ctrl_1", "Ctrl_2", "Treat_1", "Treat_2"),
#'    group = c("Control", "Control", "Treatment", "Treatment")
#' )
#' pd_obj <- create_protdata(dat = raw_data, condition = cond_df)
#'
#' # Example 1: Default mode (heatmap of individual samples)
#' p1 <- plot_proteomics_heatmap(pd_obj)
#' if (interactive()) print(p1)
#'
#' # Example 2: Summarized mode (heatmap of condition means)
#' p2 <- plot_proteomics_heatmap(pd_obj, condition = "group")
#' if (interactive()) print(p2)
#'
setGeneric("plot_proteomics_heatmap",
           function(PD, protmeta_col = NULL, genes = NULL, title = NULL, condition = NULL) {
             standardGeneric("plot_proteomics_heatmap")
           }
)

setMethod("plot_proteomics_heatmap", "ProtData",
          function(PD, protmeta_col = NULL, genes = NULL, title = NULL, condition = NULL) {
            # --- 1. Initial Data Extraction and Filtering ---
            if (!inherits(PD, "ProtData")) { stop("Error: 'PD' must be of class ProtData") }
            intensities <- PD@data
            prot_meta <- PD@prot_meta
            if (is.null(protmeta_col)) { protmeta_col <- names(prot_meta)[1] }
            if (!protmeta_col %in% colnames(prot_meta)) { stop(paste("Column", protmeta_col, "not found in PD@prot_meta")) }

            if (!is.null(genes)) {
              match_genes <- tolower(prot_meta[[protmeta_col]]) %in% tolower(genes)
              intensities <- intensities[match_genes, , drop = FALSE]
              prot_meta <- prot_meta[match_genes, , drop = FALSE]
            }
            if (nrow(intensities) == 0) { stop("No data remains after filtering.") }

            # --- 2. Prepare Data Matrix based on 'condition' argument ---
            if (!is.null(condition)) {
              # --- SUMMARIZED MODE ---
              cat("Condition provided. Summarizing replicates into means...\n")
              if (!condition %in% names(PD@condition)) { stop("'", condition, "' not found in the condition slot.") }

              # Use the summarize_replicates logic
              data_long <- intensities %>%
                tibble::rownames_to_column(var = "ProteinID") %>%
                tidyr::pivot_longer(cols = -ProteinID, names_to = "Sample", values_to = "Intensity")

              condition_to_join <- PD@condition %>% tibble::rownames_to_column(var = "Sample")

              summarized_data <- data_long %>%
                dplyr::left_join(condition_to_join, by = "Sample") %>%
                dplyr::group_by(ProteinID, .data[[condition]]) %>%
                dplyr::summarise(MeanIntensity = mean(Intensity, na.rm = TRUE)) %>%
                dplyr::ungroup()

              data_to_plot <- summarized_data %>%
                tidyr::pivot_wider(id_cols = ProteinID, names_from = .data[[condition]], values_from = MeanIntensity) %>%
                tibble::column_to_rownames(var = "ProteinID")

            } else {
              # --- INDIVIDUAL SAMPLE MODE ---
              data_to_plot <- intensities
            }

            # --- 3. Scale, Cluster, and Reshape the final data_to_plot matrix ---
            scaled_data <- t(base::scale(t(as.matrix(data_to_plot))))
            scaled_data[!is.finite(scaled_data)] <- 0
            col_clustering <- hclust(dist(t(scaled_data)))
            col_order <- colnames(scaled_data)[col_clustering$order]
            protein_labels <- make.unique(as.character(prot_meta[[protmeta_col]]))

            heatmap_data <- as.data.frame(scaled_data) %>%
              tibble::rownames_to_column(var = "Protein") %>%
              dplyr::mutate(Protein = factor(protein_labels, levels = rev(protein_labels))) %>%
              tidyr::pivot_longer(cols = -Protein, names_to = "Sample", values_to = "ZScore") %>%
              dplyr::mutate(Sample = factor(Sample, levels = col_order))

            # --- 4. Create the ggplot Heatmap ---
            plot_title <- if (!is.null(title)) title else "Proteomics Heatmap"
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
              ) + ggplot2::labs(title = plot_title)

            # Conditionally hide row labels if there are too many
            if (nrow(intensities) >= 100) {
              g <- g + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                                      axis.ticks.y = ggplot2::element_blank())
            }

            return(g)
          }
)
