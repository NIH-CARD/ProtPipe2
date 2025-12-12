#' @describeIn plot_proteomics_heatmap Method for SummarizedExperiment objects.
#' @export
# setMethod("plot_proteomics_heatmap", "SummarizedExperiment",
#           function(object, protmeta_col = NULL, genes = NULL, title = NULL, condition = NULL) {
#             # --- 1. Initial Data Extraction and Filtering ---
#             if (!inherits(object, "SummarizedExperiment")) { stop("Error: 'object' must be of class SummarizedExperiment") }
#             intensities <- assay(object) %>% as.data.frame()
#             prot_meta <- rowData(object) %>% as.data.frame()
#             if (is.null(protmeta_col)) { protmeta_col <- names(prot_meta)[1] }
#             if (!protmeta_col %in% colnames(prot_meta)) { stop(paste("Column", protmeta_col, "not found in rowData slot")) }
#
#             if (!is.null(genes)) {
#               match_genes <- tolower(prot_meta[[protmeta_col]]) %in% tolower(genes)
#               intensities <- intensities[match_genes, , drop = FALSE]
#               prot_meta <- prot_meta[match_genes, , drop = FALSE]
#             }
#             if (nrow(intensities) == 0) { stop("No data remains after filtering.") }
#
#             # --- 2. Prepare Data Matrix based on 'condition' argument ---
#             if (!is.null(condition)) {
#               # --- SUMMARIZED MODE ---
#               cat("Condition provided. Summarizing replicates into means...\n")
#               if (!condition %in% names(colData(object))) { stop("'", condition, "' not found in the condition slot.") }
#
#               # Use the summarize_replicates logic
#               data_long <- intensities %>%
#                 tibble::rownames_to_column(var = "ProteinID") %>%
#                 tidyr::pivot_longer(cols = -ProteinID, names_to = "Sample", values_to = "Intensity")
#
#               condition_to_join <- colData(object) %>%
#                 as.data.frame() %>%
#                 tibble::rownames_to_column(var = "Sample")
#
#               summarized_data <- data_long %>%
#                 dplyr::left_join(condition_to_join, by = "Sample") %>%
#                 dplyr::group_by(ProteinID, .data[[condition]]) %>%
#                 dplyr::summarise(MeanIntensity = mean(Intensity, na.rm = TRUE)) %>%
#                 dplyr::ungroup()
#
#               data_to_plot <- summarized_data %>%
#                 tidyr::pivot_wider(id_cols = ProteinID, names_from = .data[[condition]], values_from = MeanIntensity) %>%
#                 tibble::column_to_rownames(var = "ProteinID")
#
#             } else {
#               # --- INDIVIDUAL SAMPLE MODE ---
#               data_to_plot <- intensities
#             }
#
#             # --- 3. Scale, Cluster, and Reshape the final data_to_plot matrix ---
#             scaled_data <- t(base::scale(t(as.matrix(data_to_plot))))
#             scaled_data[!is.finite(scaled_data)] <- 0
#             col_clustering <- hclust(dist(t(scaled_data)))
#             col_order <- colnames(scaled_data)[col_clustering$order]
#             protein_labels <- make.unique(as.character(prot_meta[[protmeta_col]]))
#
#             heatmap_data <- as.data.frame(scaled_data) %>%
#               tibble::rownames_to_column(var = "Protein") %>%
#               dplyr::mutate(Protein = factor(protein_labels, levels = rev(protein_labels))) %>%
#               tidyr::pivot_longer(cols = -Protein, names_to = "Sample", values_to = "ZScore") %>%
#               dplyr::mutate(Sample = factor(Sample, levels = col_order))
#
#             # --- 4. Create the ggplot Heatmap ---
#             plot_title <- if (!is.null(title)) title else "Proteomics Heatmap"
#             max_val <- ceiling(max(heatmap_data$ZScore))
#             min_val <- floor(min(heatmap_data$ZScore))
#             g <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = Sample, y = Protein, fill = ZScore)) +
#               ggplot2::geom_tile(color = "black", size = 0.0) +
#               ggplot2::theme_classic() +
#               ggplot2::scale_fill_gradient2(
#                 low = "skyblue", high = "tomato1", mid = "beige",
#                 midpoint = 0,
#                 #space = "Lab", breaks = c(-max_abs_val, max_abs_val),
#                 name = "protein\nz-score",
#                 limits = c(min_val, max_val),
#               ) +
#               ggplot2::theme(
#                 axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
#                 axis.text.y = ggplot2::element_text(size = 10),
#                 axis.title = ggplot2::element_blank()
#               ) + ggplot2::labs(title = plot_title)
#
#             # Conditionally hide row labels if there are too many
#             if (nrow(intensities) >= 100) {
#               g <- g + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
#                                       axis.ticks.y = ggplot2::element_blank())
#             }
#
#             return(g)
#           }
# )
setMethod("plot_proteomics_heatmap", "SummarizedExperiment",
          function(object, protmeta_col = NULL, genes = NULL, title = NULL, condition = NULL,
                   cluster_cols = FALSE, cluster_rows = FALSE) { # <--- Added parameters here

            # --- 1. Initial Data Extraction and Filtering ---
            if (!inherits(object, "SummarizedExperiment")) { stop("Error: 'object' must be of class SummarizedExperiment") }
            intensities <- assay(object) %>% as.data.frame()
            prot_meta <- rowData(object) %>% as.data.frame()
            if (is.null(protmeta_col)) { protmeta_col <- names(prot_meta)[1] }
            if (!protmeta_col %in% colnames(prot_meta)) { stop(paste("Column", protmeta_col, "not found in rowData slot")) }

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
              if (!condition %in% names(colData(object))) { stop("'", condition, "' not found in the condition slot.") }

              # Use the summarize_replicates logic
              data_long <- intensities %>%
                tibble::rownames_to_column(var = "ProteinID") %>%
                tidyr::pivot_longer(cols = -ProteinID, names_to = "Sample", values_to = "Intensity")

              condition_to_join <- colData(object) %>%
                as.data.frame() %>%
                tibble::rownames_to_column(var = "Sample")

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

            # -- Conditional Column Clustering --
            if (cluster_cols) {
              col_clustering <- hclust(dist(t(scaled_data)))
              col_order <- colnames(scaled_data)[col_clustering$order]
            } else {
              col_order <- colnames(scaled_data)
            }

            protein_labels <- make.unique(as.character(prot_meta[[protmeta_col]]))

            # -- Conditional Row Clustering --
            if (cluster_rows) {
              row_clustering <- hclust(dist(scaled_data))
              # Order based on dendrogram
              row_levels <- protein_labels[row_clustering$order]
            } else {
              # Default: Reverse order so the first row in matrix appears at the top of the plot
              row_levels <- rev(protein_labels)
            }

            heatmap_data <- as.data.frame(scaled_data) %>%
              tibble::rownames_to_column(var = "Protein") %>%
              # Apply the conditional row levels
              dplyr::mutate(Protein = factor(protein_labels, levels = row_levels)) %>%
              tidyr::pivot_longer(cols = -Protein, names_to = "Sample", values_to = "ZScore") %>%
              # Apply the conditional column levels
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
