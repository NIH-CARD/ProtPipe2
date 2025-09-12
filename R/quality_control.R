#' Get Protein Group Counts Per Sample
#'
#' @description
#' Calculates the number of non-missing protein groups (or features) for each
#' sample in a ProtData object.
#'
#' @param PD A `ProtData` object.
#'
#' @return A data frame with two columns: `Sample` (the sample names) and `N`
#'   (the count of non-NA proteins for that sample).
#'
#' @export
#'
#' @examples
#' # Create a sample ProtData object with some missing values
#' raw_data <- data.frame(
#'   Gene = c("GENEA", "GENEB", "GENEC", "GENED"),
#'   SampleA = c(10, 11, NA, 13),
#'   SampleB = c(12, 13, 14, 15),
#'   SampleC = c(NA, NA, 13, 14)
#' )
#' pd_obj <- create_protdata(dat = raw_data)
#'
#' # Get the protein counts for each sample
#' counts <- get_pg_counts(pd_obj)
#' print(counts)
#' # Expected output:
#' #   Sample N
#' # 1 SampleA 3
#' # 2 SampleB 4
#' # 3 SampleC 2
#'
setGeneric("get_pg_counts", function(PD) standardGeneric("get_pg_counts"))

setMethod("get_pg_counts", "ProtData",
          function(PD){
            # get the number of protein groups per sample
            N_values <- colSums(!is.na(getData(PD)))
            pgcounts <- data.frame(Sample = names(N_values), N = N_values)

            return(pgcounts)
          }
)


#' Plot Protein Group Counts
#'
#' @description
#' Generates a bar plot showing the number of identified protein groups per sample.
#' The plot can either show counts for each individual sample, or show the
#' summarized mean and standard deviation for samples grouped by a condition.
#'
#' @param PD A `ProtData` object.
#' @param condition Optional. A character string specifying a column name in the
#'   `condition` slot. If provided, samples will be grouped by this variable and
#'   the plot will show the mean and standard deviation of counts per group.
#'   If `NULL` (the default), each sample is plotted individually.
#'
#' @return A `ggplot` object.
#'
#' @export
#'
#' @examples
#' # Create sample data
#' raw_data <- data.frame(
#'   Gene = c("GENEA", "GENEB", "GENEC", "GENED"),
#'   Ctrl_1 = c(10, 11, NA, 13),
#'   Ctrl_2 = c(12, 13, 14, 15),
#'   Treat_1 = c(NA, NA, 13, 14),
#'   Treat_2 = c(10, 11, 12, 13)
#' )
#' cond_df <- data.frame(
#'    SampleID = c("Ctrl_1", "Ctrl_2", "Treat_1", "Treat_2"),
#'    group = c("Control", "Control", "Treatment", "Treatment")
#' )
#' pd_obj <- create_protdata(dat = raw_data, condition = cond_df)
#'
#' # Example 1: Plot individual sample counts
#' p1 <- plot_pg_counts(pd_obj)
#' if (interactive()) {
#'   print(p1)
#' }
#'
#' # Example 2: Plot counts grouped by the 'group' condition
#' p2 <- plot_pg_counts(pd_obj, condition = "group")
#' if (interactive()) {
#'   print(p2)
#' }
#'
setGeneric("plot_pg_counts", function(PD, condition = NULL) standardGeneric("plot_pg_counts"))

setMethod("plot_pg_counts", "ProtData",
          function(PD, condition = NULL) {

            pgcounts <- get_pg_counts(PD)

            # Order samples by ascending counts
            n_samples <- nrow(pgcounts)
            if (is.null(condition)){
              if (n_samples > 20) {
                p=ggplot2::ggplot(pgcounts, ggplot2::aes(x=Sample, y=N)) +
                  ggplot2::geom_bar(stat="identity", fill="#67a9cf")+
                  ggplot2::theme_classic()+
                  ggplot2::labs(fill = "",x="",y='Number of Protein Groups')+
                  ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90))
              }
              if (n_samples <= 20) {
                p=ggplot2::ggplot(pgcounts, ggplot2::aes(x=Sample, y=N)) +
                  ggplot2::geom_bar(stat="identity", fill="#67a9cf")+
                  ggplot2::theme_classic()+
                  ggplot2::labs(fill = "",x="",y='Number of Protein Groups')+
                  ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90))+
                  ggplot2::geom_text(ggplot2::aes(label=N, y=N + (0.05*max(pgcounts$N))))
              }
            }else{
              # group by condition
              condition_file <- PD@condition
              if (condition %in% colnames(condition_file)){
                condition_file <- condition_file %>%
                  dplyr::mutate(Sample = rownames(condition_file)) %>%
                  dplyr::select(c(!!rlang::sym(condition), Sample))
                pgcounts <- pgcounts %>%
                  dplyr::left_join(condition_file, by = "Sample") %>%
                  dplyr::rename(Condition = !!rlang::sym(condition))  # Rename the dynamic column to 'Condition'
                summary_data <- pgcounts %>%
                  dplyr::group_by(Condition) %>%
                  dplyr::summarize(mean = mean(N), sd = sd(N)) %>%
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

#' Plot Boxplots of Sample Intensity Distributions
#'
#' @description
#' Generates boxplots to visualize and compare the distribution of protein
#' intensities for each sample. This is a common quality control plot to check
#' for normalization issues or identify potential outlier samples.
#'
#' @param PD A `ProtData` object.
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
#' pd_obj <- create_protdata(dat = raw_data)
#'
#' # Generate the boxplot of intensities
#' p <- plot_pg_intensities(pd_obj)
#' if (interactive()) {
#'   print(p)
#' }
#'
setGeneric("plot_pg_intensities", function(PD) standardGeneric("plot_pg_intensities"))

setMethod("plot_pg_intensities", "ProtData",
          function(PD) {
            # Assuming PD@data is a data frame with proteins as rows and samples as columns
            dat <- PD@data  # Or however you access your data frame in the wide format
            dat_long <- reshape2::melt(dat,
                                       measure.vars=names(dat)[sapply(dat, function(x) all(is.numeric(x)))],
                                       variable.name='Sample',
                                       value.name='Intensity')
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

#' Calculate Coefficient of Variation (CV) for Protein Groups
#'
#' @description
#' Calculates the coefficient of variation (CV) for each protein, grouped by an
#' experimental condition. The CV is a standardized measure of dispersion,
#' calculated as the standard deviation divided by the mean, and is often used
#' to assess the reproducibility of replicates within a condition group.
#'
#' @param PD A `ProtData` object.
#' @param condition A character string specifying the column name in the
#'   `condition` slot which defines the replicate groups.
#' @param min_samples The minimum number of replicates in a group required to
#'   calculate CVs. Groups with fewer samples will be skipped. Defaults to 2.
#'
#' @return A data frame in long format with columns for `Protein`, `CV`
#'   (in percent), and `Condition`.
#'
#' @export
#'
#' @examples
#' # Create sample data with two conditions, each with 3 replicates
#' raw_data <- data.frame(
#'   Protein.ID = paste0("Prot", 1:5),
#'   Gene = paste0("GENE", LETTERS[1:5]),
#'   Ctrl_1 = rnorm(5, 10, 1), Ctrl_2 = rnorm(5, 10, 1.5), Ctrl_3 = rnorm(5, 10, 1.2),
#'   Treat_1 = rnorm(5, 15, 2), Treat_2 = rnorm(5, 15, 2.5), Treat_3 = rnorm(5, 15, 2.2)
#' )
#' cond_df <- data.frame(
#'    SampleID = c("Ctrl_1", "Ctrl_2", "Ctrl_3", "Treat_1", "Treat_2", "Treat_3"),
#'    group = c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment")
#' )
#' pd_obj <- create_protdata(dat = raw_data, condition = cond_df)
#'
#' # Calculate the CVs for each protein within the 'group' condition
#' cv_results <- get_CVs(pd_obj, condition = "group")
#' head(cv_results)
#'
setGeneric("get_CVs", function(PD, condition, min_samples = 2) standardGeneric("get_CVs"))

setMethod("get_CVs", "ProtData",
          function(PD, condition, min_samples = 2) {
            condition_file <- PD@condition
            if (!condition %in% colnames(condition_file)) {
              stop("The selected condition does not appear in the condition file.")
            }

            intensities <- as.matrix(PD@data)
            if (!is.numeric(intensities)) {
              stop("PD@data must contain only numeric values.")
            }

            conds <- condition_file[[condition]]
            unique_conds <- (unique(conds))

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
                Protein = PD@prot_meta[[1]],
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

#' Plot Coefficient of Variation (CV) Distributions
#'
#' @description
#' Generates a violin or jitter plot to visualize the distribution of protein
#' CVs for each experimental condition. This plot is useful for comparing the
#' reproducibility and technical variance across different sample groups.
#'
#' @param PD A `ProtData` object.
#' @param condition A character string specifying the column name in the
#'   `condition` slot which defines the replicate groups.
#' @param plot_type The type of plot to generate. Must be either `"violin"`
#'   (the default, which also includes a narrow boxplot) or `"jitter"`.
#'
#' @return A `ggplot` object.
#'
#' @export
#'
#' @examples
#' # Create sample data with two conditions, each with 3 replicates
#' raw_data <- data.frame(
#'   Protein.ID = paste0("Prot", 1:10),
#'   Ctrl_1 = rnorm(10, 10, 1), Ctrl_2 = rnorm(10, 10, 1.5), Ctrl_3 = rnorm(10, 10, 1.2),
#'   Treat_1 = rnorm(10, 15, 3), Treat_2 = rnorm(10, 15, 3.5), Treat_3 = rnorm(10, 15, 3.2)
#' )
#' cond_df <- data.frame(
#'    SampleID = c("Ctrl_1", "Ctrl_2", "Ctrl_3", "Treat_1", "Treat_2", "Treat_3"),
#'    group = c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment")
#' )
#' pd_obj <- create_protdata(dat = raw_data, condition = cond_df)
#'
#' # Generate the default violin plot
#' p_violin <- plot_CVs(pd_obj, condition = "group")
#' if (interactive()) {
#'   print(p_violin)
#' }
#'
#' # Generate the jitter plot instead
#' p_jitter <- plot_CVs(pd_obj, condition = "group", plot_type = "jitter")
#' if (interactive()) {
#'   print(p_jitter)
#' }
#'
setGeneric("plot_CVs",
           function(PD, condition, plot_type = "violin") {
             standardGeneric("plot_CVs")
           }
)

setMethod("plot_CVs", "ProtData",
          function(PD, condition, plot_type = "violin") {
            plot_type <- match.arg(plot_type, choices = c("violin", "jitter"))

            cv_df <- get_CVs(PD, condition)
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

#' Calculate Pairwise Sample Correlations
#'
#' @description
#' Computes a pairwise correlation matrix for all samples based on their protein
#' intensity profiles. The result is returned in a long "tidy" format data table.
#'
#' @details
#' The function name is `get_spearman`, but it can compute other correlation
#' coefficients by changing the `method` argument.
#'
#' @param PD A `ProtData` object.
#' @param method The correlation method to use. Can be one of `"spearman"` (the
#'   default), `"pearson"`, or `"kendall"`.
#'
#' @return A `data.table` with three columns: `SampleA`, `SampleB`, and a third
#'   column named for the method used (e.g., `Spearman`), containing the
#'   pairwise correlation coefficients.
#'
#' @export
#'
#' @examples
#' # Create sample data
#' raw_data <- data.frame(
#'   Gene = c("GENEA", "GENEB", "GENEC", "GENED"),
#'   SampleA = c(10, 20, 15, 12),
#'   SampleB = c(11, 21, 16, 13), # Correlated with A
#'   SampleC = c(25, 10, 30, 5)   # Not correlated with A/B
#' )
#' pd_obj <- create_protdata(dat = raw_data)
#'
#' # Get the spearman correlation table
#' corr_table <- get_spearman(pd_obj)
#' print(corr_table)
#'
setGeneric("get_spearman", function(PD, method = 'spearman') standardGeneric("get_spearman"))

setMethod("get_spearman", "ProtData",
          function(PD, method = 'spearman') {
            DT <- getData(PD)
            #### Pairwise correlations between sample columns
            #dt.samples <- DT[,-c(1:2)]      # Ignore info columns (subset to only intensity values)
            dt.samples <- DT[, sapply(DT, is.numeric)] #better way of getting just intensity columns

            dt.corrs <- cor(as.matrix(na.omit(dt.samples)+1), method=method)

            # Format correlations as 3 digits
            dt.corrs <- data.table::data.table(reshape2::melt(dt.corrs, measure.vars=dt.corrs[,rn], value.name='Spearman'))
            dt.corrs <- dt.corrs[! is.na('Spearman')]
            data.table::setnames(dt.corrs, c('Var1', 'Var2'), c('SampleA','SampleB'))
            dt.corrs <- dt.corrs %>% dplyr::mutate(Spearman = round(Spearman, 3))

            return(dt.corrs[])
          }
)

#' Plot a Sample Correlation Heatmap
#'
#' @description
#' Creates a heatmap of pairwise sample correlations, calculated via the Spearman
#' method. The function includes a "smart sorting" feature to correctly order
#' axes based on numeric parts of labels (e.g., "Day1", "Day2", "Day10").
#'
#' @param PD A `ProtData` object.
#' @param order_by Optional. A character string specifying a column in the
#'   condition slot by which to order the samples on the heatmap axes. If the
#'   column contains numeric parts, a natural sort is applied.
#' @param label_by Optional. A character string specifying a column in the
#'   condition slot to use for relabeling the heatmap axes.
#'
#' @return A `ggplot` object representing the heatmap.
#'
#' @export
#'
#' @examples
#' # Create sample data with timepoints for sorting
#' raw_data <- data.frame(
#'   Gene = paste0("G", 1:5),
#'   Sample_D10 = rnorm(5, 15), # Treatment
#'   Sample_D1 = rnorm(5, 10),  # Control
#'   Sample_D2 = rnorm(5, 10)   # Control
#' )
#' cond_df <- data.frame(
#'   SampleID = c("Sample_D1", "Sample_D2", "Sample_D10"),
#'   Timepoint = c("Day1", "Day2", "Day10"),
#'   Group = c("Control", "Control", "Treatment")
#' )
#' pd_obj <- create_protdata(dat = raw_data, condition = cond_df)
#'
#' # Example 1: Default alphanumeric sorting (D1, D10, D2)
#' p1 <- plot_correlation_heatmap(pd_obj)
#' if (interactive()) print(p1)
#'
#' # Example 2: Smart sorting by 'Timepoint' (D1, D2, D10)
#' p2 <- plot_correlation_heatmap(pd_obj, order_by = "Timepoint")
#' if (interactive()) print(p2)
#'
#' # Example 3: Smart sorting and relabeling by 'Group'
#' p3 <- plot_correlation_heatmap(pd_obj, order_by = "Timepoint", label_by = "Group")
#' if (interactive()) print(p3)
#'
setGeneric("plot_correlation_heatmap",
           function(PD, order_by = NULL, label_by = NULL) {
             standardGeneric("plot_correlation_heatmap")
           }
)

setMethod("plot_correlation_heatmap", "ProtData",
          function(PD, order_by = NULL, label_by = NULL) {

            # --- 1. Calculate Correlation Data ---
            DT.corrs <- get_spearman(PD)

            # --- 2. Determine Sample Order ---
            metadata <- getCondition(PD)

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
              ggplot2::geom_tile(color = "white", size = 0.5) +
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
