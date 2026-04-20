#' @rdname get_pg_counts
#' @export
setMethod("get_pg_counts", "SummarizedExperiment",
          function(object){
            # get the number of protein groups per sample
            n_values <- colSums(!is.na(assay(object)))
            pgcounts <- data.frame(Sample = names(n_values), Protein_Groups = n_values)

            return(pgcounts)
          }
)

#' @rdname plot_pg_counts
#' @export
setMethod("plot_pg_counts", "SummarizedExperiment",
          function(object, condition = NULL) {

            pgcounts <- get_pg_counts(object)

            # Order samples by ascending counts
            n_samples <- nrow(pgcounts)
            if (is.null(condition)){
              if (n_samples > 20) {
                p=ggplot2::ggplot(pgcounts, ggplot2::aes(x=Sample, y=Protein_Groups)) +
                  ggplot2::geom_bar(stat="identity", fill="#67a9cf")+
                  ggplot2::theme_classic()+
                  ggplot2::labs(fill = "",x="",y='Number of Protein Groups')+
                  ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90))
              }
              if (n_samples <= 20) {
                p=ggplot2::ggplot(pgcounts, ggplot2::aes(x=Sample, y=Protein_Groups)) +
                  ggplot2::geom_bar(stat="identity", fill="#67a9cf")+
                  ggplot2::theme_classic()+
                  ggplot2::labs(fill = "",x="",y='Number of Protein Groups')+
                  ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90))+
                  ggplot2::geom_text(ggplot2::aes(label=Protein_Groups, y=Protein_Groups + (0.05*max(pgcounts$Protein_Groups))))
              }
            }else{
              # group by condition
              condition_file <- as.data.frame(colData(object))
              if (condition %in% colnames(condition_file)){
                condition_file <- condition_file %>%
                  dplyr::mutate(Sample = rownames(condition_file)) %>%
                  dplyr::select(c(!!rlang::sym(condition), Sample))
                pgcounts <- pgcounts %>%
                  dplyr::left_join(condition_file, by = "Sample") %>%
                  dplyr::rename(Condition = !!rlang::sym(condition))  # Rename the dynamic column to 'Condition'
                summary_data <- pgcounts %>%
                  dplyr::group_by(Condition) %>%
                  dplyr::summarize(mean = mean(Protein_Groups), sd = sd(Protein_Groups)) %>%
                  dplyr::arrange(Condition)
                p=ggplot2::ggplot(summary_data, ggplot2::aes(x=as.factor(Condition), y=mean)) +
                  ggplot2::geom_bar(stat="identity",fill="#67a9cf", position= ggplot2::position_dodge())+
                  ggplot2::theme_classic()+
                  ggplot2::geom_errorbar(ggplot2::aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                                         position=ggplot2::position_dodge(.9))+
                  ggplot2::labs(fill = "",x="",y='Number of Protein Groups')
              }else{
                stop("the selected condition does not appear in the condition file")
              }
            }
            return(p)
          }
)

#' @rdname plot_pg_intensities
#' @export
setMethod("plot_pg_intensities", "SummarizedExperiment",
          function(object) {
            # Assuming PD@data is a data frame with proteins as rows and samples as columns
            dat <- as.data.frame(assay(object))  # Or however you access your data frame in the wide format
            dat_long <- reshape2::melt(dat,
                                       measure.vars=names(dat)[sapply(dat, function(x) all(is.numeric(x)))],
                                       variable.name='Sample',
                                       value.name='Intensity')
            dat_long <<-dat_long
            dat_long <- dat_long[dat_long$Intensity>0,]
            dat_long <- dat_long[rowSums(is.na(dat_long)) < ncol(dat_long),]
            # Plot the boxplot
            g <- ggplot2::ggplot(dat_long, ggplot2::aes(x = Sample, y = log10(Intensity))) +
              ggplot2::geom_boxplot(outlier.shape = NA, fill = "#67a9cf") +
              ggplot2::theme_classic() +
              ggplot2::labs(fill = "", x = "", y = "Log10 Protein Intensity") +
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
              ggplot2::geom_boxplot(width = 0.1) +
              ggplot2::geom_hline(color = '#ef8a62', linetype = 'dashed',
                                  ggplot2::aes(yintercept = quantile(log10(Intensity), 0.50, na.rm = TRUE)))

            return(g)
          }
)

#' @rdname get_CVs
#' @export
setMethod("get_CVs", "SummarizedExperiment",
          function(object, condition, min_samples = 2) {
            condition_file <- as.data.frame(colData(object))
            if (!condition %in% colnames(condition_file)) {
              stop("The selected condition does not appear in the condition file.")
            }

            intensities <- as.matrix(assay(object))
            if (!is.numeric(intensities)) {
              stop("The assay must contain only numeric values.")
            }

            conds <- condition_file[[condition]]
            unique_conds <- (unique(conds))
            group_sizes <- table(conds, useNA = "no")
            undersized_groups <- names(group_sizes[group_sizes < min_samples])

            if (length(undersized_groups) > 0) {
              stop(
                paste0(
                  "Selected groups must each contain at least ", min_samples,
                  " samples. Too few samples in: ",
                  paste(undersized_groups, collapse = ", "), "."
                )
              )
            }

            cv_list <- lapply(unique_conds, function(cond) {
              idx <- which(conds == cond)
              if (length(idx) < min_samples) {
                return(NULL)
              }

              sub_data <- intensities[, idx, drop = FALSE]
              means <- matrixStats::rowMeans2(sub_data, na.rm = TRUE)
              sds <- matrixStats::rowSds(sub_data, na.rm = TRUE)
              cvs <- sds * 100/ means

              data.frame(
                Protein = rowData(object)[[1]],
                CV = cvs,
                Condition = cond,
                stringsAsFactors = FALSE
              )
            })

            cv_df <- do.call(rbind, cv_list)
            cv_df <- cv_df[!is.na(cv_df$CV), ]
            return(cv_df)
          }
)

#' @rdname plot_CVs
#' @export
setMethod("plot_CVs", "SummarizedExperiment",
          function(object, condition, plot_type = "violin") {
            plot_type <- match.arg(plot_type, choices = c("violin", "jitter"))

            cv_df <- get_CVs(object, condition)
            cv_df$Condition <- as.factor(cv_df$Condition)
            p <- ggplot2::ggplot(cv_df, ggplot2::aes(x = Condition, y = CV)) +
              ggplot2::theme_classic() +
              ggplot2::labs(x = "", y = "Coefficient of Variation (%)") +
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
              ggplot2::geom_hline(
                ggplot2::aes(yintercept = stats::quantile(CV, 0.5, na.rm = TRUE)),
                color = "#ef8a62", linetype = "dashed"
              )

            if (plot_type == "violin") {
              p <- p +
                ggplot2::geom_violin(fill = "#67a9cf", color = "black", trim = FALSE) +
                ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA)
            } else if (plot_type == "jitter") {
              p <- p +
                ggplot2::geom_jitter(width = 0.15, alpha = 0.5, size = 0.7)
            }

            return(p)
          }
)

#' @rdname get_sample_correlation
#' @export
setMethod("get_sample_correlation", "SummarizedExperiment",
          function(object, method = 'spearman', num_features = NULL) {
            dt.samples <- as.matrix(assay(object))

            if (!is.numeric(dt.samples)) {
              stop("The assay must contain only numeric values.")
            }

            if (!is.null(num_features)) {
              if (!is.numeric(num_features) || length(num_features) != 1 || is.na(num_features) || num_features < 1) {
                stop("num_features must be a single positive number.")
              }

              n_features <- min(as.integer(num_features), nrow(dt.samples))
              row_variances <- matrixStats::rowVars(dt.samples, na.rm = TRUE)
              row_variances[is.na(row_variances)] <- -Inf
              selected_rows <- order(row_variances, decreasing = TRUE)[seq_len(n_features)]
              dt.samples <- dt.samples[selected_rows, , drop = FALSE]
            }

            dt.corrs <- cor(dt.samples + 1,
                            method = method,
                            use = "pairwise.complete.obs")

            # Format correlations as 3 digits
            dt.corrs <- data.table::data.table(reshape2::melt(dt.corrs, measure.vars=dt.corrs[,rn], value.name='Spearman'))
            dt.corrs <- dt.corrs[! is.na('Spearman')]
            data.table::setnames(dt.corrs, c('Var1', 'Var2'), c('SampleA','SampleB'))
            dt.corrs <- dt.corrs %>% dplyr::mutate(Spearman = round(Spearman, 3))

            return(dt.corrs[])
          }
)

#' @rdname plot_correlation_heatmap
#' @export
setMethod("plot_correlation_heatmap", "SummarizedExperiment",
          function(object, order_by = NULL, label_by = NULL, num_features = NULL) {

            # --- 1. Calculate Correlation Data ---
            DT.corrs <- get_sample_correlation(object, num_features = num_features)

            # --- 2. Determine Sample Order ---
            metadata <- as.data.frame(colData(object))

            # Default order is sorted sample names
            sample_order <- sort(rownames(metadata))

            if (!is.null(order_by)) {
              if (!order_by %in% names(metadata)) {
                stop("'", order_by, "' is not a valid column in the condition data.")
              }

              # -------------------------------------------------------------------- #
              # --- NEW: Smart Sorting Logic ---
              # -------------------------------------------------------------------- #
              ordering_vector <- metadata[[order_by]]

              # Attempt to extract all digits from the ordering variable
              # The suppressWarnings (@) is used in case there are no numbers, to avoid a warning message.
              numeric_part <- suppressWarnings(as.numeric(gsub("\\D", "", ordering_vector)))

              # If the numeric conversion was successful (i.e., we found numbers), sort by them.
              if (!all(is.na(numeric_part))) {
                cat("Numeric values detected in ordering variable. Applying natural sort.\n")
                # Order the sample names based on the extracted numeric values
                sample_order <- rownames(metadata)[order(numeric_part)]
              } else {
                # Otherwise, if no numbers were found, fall back to the default alphanumeric sort.
                cat("No numeric values detected in ordering variable. Applying alphanumeric sort.\n")
                sample_order <- rownames(metadata)[order(ordering_vector)]
              }
              # -------------------------------------------------------------------- #
            }

            # --- 3. Apply the Order to the Data ---
            DT.corrs$SampleA <- factor(DT.corrs$SampleA, levels = sample_order)
            DT.corrs$SampleB <- factor(DT.corrs$SampleB, levels = sample_order)

            # --- 4. Define Plot Aesthetics ---
            max_limit <- max(DT.corrs$Spearman)
            min_limit <- min(DT.corrs$Spearman)
            mid_limit <- round((max_limit + min_limit) / 2, 3)

            # --- 5. Create the Plot (Single, clean ggplot call) ---
            g <- ggplot2::ggplot(DT.corrs, ggplot2::aes(x = SampleA, y = SampleB, fill = Spearman)) +
              ggplot2::geom_tile(color = "white", linewidth = 0.5) +
              ggplot2::theme_classic() +
              ggplot2::scale_fill_gradient2(
                low = "skyblue", high = "tomato1", mid = "white",
                midpoint = mid_limit, limit = c(min_limit, max_limit),
                space = "Lab", breaks = c(min_limit, mid_limit, max_limit),
                name = "Spearman\nCorrelation"
              ) +
              ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
                axis.text.y = ggplot2::element_text(size = 10),
                axis.title = ggplot2::element_blank()
              )

            # --- 6. Conditionally Add Axis Labels ---
            if (!is.null(label_by)) {
              if (!label_by %in% names(metadata)) {
                stop("'", label_by, "' is not a valid column in the condition data.")
              }
              label_map <- setNames(as.character(metadata[[label_by]]), rownames(metadata))

              g <- g +
                ggplot2::scale_x_discrete(labels = label_map) +
                ggplot2::scale_y_discrete(labels = label_map)
            }

            return(g)
          }
)
