#' # --- 1. Required Packages ---
#' # Make sure you have these packages installed for the function to work
#' # install.packages(c("dplyr", "ggplot2", "ggpubr", "rlang"))
#'
#' # --- 2. Generic Function Definition ---
#' #' Generate a bar chart of protein intensity values
#' #'
#' #' Creates a ggplot bar chart comparing the intensity of a single protein
#' #' either across all samples or grouped by a condition.
#' #'
#' #' @param object The `SummarizedExperiment` object.
#' #' @param prot A string specifying the name of the protein to plot.
#' #' @param prot_meta_col A string naming the column in the `rowData` slot to search for the protein. Defaults to the first column.
#' #' @param condition (Optional) A string naming the column in the `colData` slot to group samples by.
#' #' @return A ggplot object.
#' #' @export
#' setGeneric("compare_protein",
#'            def = function(object, prot, prot_meta_col = NULL, condition = NULL) {
#'              standardGeneric("compare_protein")
#'            }
#' )
#'
#' # --- 3. S4 Method Implementation ---
#' #' @describeIn compare_protein Method for ProtData objects.
#' setMethod("compare_protein",
#'           signature(object = "SummarizedExperiment"),
#'           function(object, prot, prot_meta_col = NULL, condition = NULL) {
#'
#'             # --- Input Validation and Data Extraction ---
#'
#'             if (length(prot) != 1 || !is.character(prot)) {
#'               stop("'prot' must be a single character string.")
#'             }
#'
#'             # 1. Validate and default the prot_meta_col
#'             if (is.null(prot_meta_col)) {
#'               prot_meta_col <- names(rowData(object))[[1]]
#'               message("`prot_meta_col` not provided. Using the first column: '", prot_meta_col, "'")
#'             }
#'             if (!prot_meta_col %in% names(rowData(object))) {
#'               stop("'", prot_meta_col, "' not found in the prot_meta slot of the object.")
#'             }
#'
#'             # 2. Find the protein row index
#'             # Use case-insensitive matching
#'             match_indices <- which(tolower(rowData(object)[[prot_meta_col]]) == tolower(prot))
#'
#'             if (length(match_indices) == 0) {
#'               stop("Analyte '", prot, "' not found in column '", prot_meta_col, "'.")
#'             }
#'
#'             # 3. Handle multiple matches
#'             if (length(match_indices) > 1) {
#'               warning("Multiple matches found for '", prot, "'. Using the first one.")
#'               match_indices <- match_indices[1]
#'             }
#'
#'             # 4. Extract the intensity data for the single protein
#'             intensity_vector <- assay(object)[match_indices, ]
#'
#'             # 5. Create a tidy data frame for plotting
#'             plotting_df <- data.frame(
#'               SampleID = names(intensity_vector),
#'               Intensity = as.numeric(intensity_vector)
#'             )
#'
#'
#'             # --- Plotting Logic ---
#'
#'             # Case 1: No condition provided, plot individual samples
#'             if (is.null(condition)) {
#'
#'               g <- ggplot2::ggplot(plotting_df, ggplot2::aes(x = SampleID, y = Intensity)) +
#'                 ggplot2::geom_col(fill = "steelblue", color = "black") +
#'                 ggplot2::labs(
#'                   title = paste("Intensity for Protein:", prot),
#'                   x = "Sample",
#'                   y = "Intensity"
#'                 ) +
#'                 ggplot2::theme_classic() +
#'                 ggplot2::theme(
#'                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
#'                 )
#'
#'               # Case 2: Condition provided, plot grouped means with stats
#'             }else{
#'               if (!condition %in% names(colData(object))) { stop("'", condition, "' not found in colData slot.") }
#'
#'               condition_df <- colData(object) %>%
#'                 as.data.frame() %>%
#'                 dplyr::mutate(SampleID = as.character(rownames(.)))
#'
#'               plotting_df <- dplyr::left_join(plotting_df, condition_df, by = "SampleID")
#'
#'               groups <- unique(na.omit(plotting_df[[condition]]))
#'
#'               # Main plot call
#'               g <- ggplot2::ggplot(plotting_df,
#'                                    ggplot2::aes(x = !!rlang::sym(condition), y = Intensity, fill = !!rlang::sym(condition))
#'               ) +
#'                 # Draw the bars FIRST
#'                 ggplot2::stat_summary(fun = mean, geom = "bar", color = "black") +
#'                 # Draw ALL error bars SECOND, so they are on top.
#'                 ggplot2::stat_summary(
#'                   fun.data = ggplot2::mean_sdl, fun.args = list(mult = 1),
#'                   geom = "errorbar", width = 0.2, color = "black"
#'                 ) +
#'                 ggplot2::labs(
#'                   title = paste("Mean Intensity for Protein:", prot),
#'                   subtitle = paste("Grouped by:", condition),
#'                   x = condition,
#'                   y = "Mean Intensity (+/- SD)"
#'                 ) +
#'                 ggplot2::theme_classic() +
#'                 ggplot2::theme(legend.position = "none")
#'
#'               # Pre-calculate stats and filter for ONLY significant pairs
#'               if (length(groups) >= 2) {
#'                 stat_test <- tryCatch({
#'
#'                   # --- This is the code that might fail ---
#'                   ggpubr::compare_means(as.formula(paste("Intensity ~", condition)), data = plotting_df) %>%
#'                     dplyr::filter(p.signif != "ns")
#'                   # ----------------------------------------
#'
#'                 }, error = function(e) {
#'
#'                   # If an error occurs, stop the execution and provide a clear message.
#'                   # We also print the original error ('e$message') for easier debugging if needed.
#'                   stop(paste("Statistical comparison failed. This can happen if a group has too few non-missing observations to perform the test.\nOriginal error:", e$message))
#'
#'                 })
#'
#'                 # Only add the layer if there are significant results to show
#'                 if (nrow(stat_test) > 0) {
#'                   # Convert the significant pairs data.frame into the list format ggpubr requires
#'                   significant_comparisons <- lapply(1:nrow(stat_test), function(i) {
#'                     c(stat_test$group1[i], stat_test$group2[i])
#'                   })
#'
#'                   g <- g + ggpubr::stat_compare_means(
#'                     comparisons = significant_comparisons,
#'                     label = "p.signif"
#'                   )
#'                 }
#'               }
#'             }
#'
#'             return(g)
#'           })
#'
# --- 2. Generic Function Definition ---
#' Generate a bar chart of protein intensity values
#'
#' Creates a ggplot bar chart comparing the intensity of a single protein
#' either across all samples or grouped by a condition.
#'
#' @param object The `SummarizedExperiment` object.
#' @param prot A string specifying the name of the protein to plot.
#' @param prot_meta_col A string naming the column in the `rowData` slot to search for the protein. Defaults to the first column.
#' @param condition (Optional) A string naming the column in the `colData` slot to group samples by.
#' @param selected_groups (Optional) A character vector of group names to filter by. Only used if 'condition' is provided.
#' @return A ggplot object.
#' @export
setGeneric("compare_protein",
           def = function(object, prot, prot_meta_col = NULL, condition = NULL, selected_groups = NULL) {
             standardGeneric("compare_protein")
           }
)

# --- 3. S4 Method Implementation ---
#' @describeIn compare_protein Method for SummarizedExperiment objects.
setMethod("compare_protein",
          signature(object = "SummarizedExperiment"),
          function(object, prot, prot_meta_col = NULL, condition = NULL, selected_groups = NULL) {

            # --- Input Validation and Data Extraction ---

            if (length(prot) != 1 || !is.character(prot)) {
              stop("'prot' must be a single character string.")
            }

            # 1. Validate and default the prot_meta_col
            if (is.null(prot_meta_col)) {
              prot_meta_col <- names(rowData(object))[[1]]
              message("`prot_meta_col` not provided. Using the first column: '", prot_meta_col, "'")
            }
            if (!prot_meta_col %in% names(rowData(object))) {
              stop("'", prot_meta_col, "' not found in the prot_meta slot of the object.")
            }

            # 2. Find the protein row index
            # Use case-insensitive matching
            match_indices <- which(tolower(rowData(object)[[prot_meta_col]]) == tolower(prot))

            if (length(match_indices) == 0) {
              stop("Analyte '", prot, "' not found in column '", prot_meta_col, "'.")
            }

            # 3. Handle multiple matches
            if (length(match_indices) > 1) {
              warning("Multiple matches found for '", prot, "'. Using the first one.")
              match_indices <- match_indices[1]
            }

            # 4. Extract the intensity data for the single protein
            intensity_vector <- assay(object)[match_indices, ]

            # 5. Create a tidy data frame for plotting
            plotting_df <- data.frame(
              SampleID = names(intensity_vector),
              Intensity = as.numeric(intensity_vector)
            )


            # --- Plotting Logic ---

            # Case 1: No condition provided, plot individual samples
            if (is.null(condition)) {

              # Warn if selected_groups provided but no condition
              if(!is.null(selected_groups)) {
                warning("'selected_groups' argument ignored because 'condition' is NULL.")
              }

              g <- ggplot2::ggplot(plotting_df, ggplot2::aes(x = SampleID, y = Intensity)) +
                ggplot2::geom_col(fill = "steelblue", color = "black") +
                ggplot2::labs(
                  title = paste("Intensity for Protein:", prot),
                  x = "Sample",
                  y = "Intensity"
                ) +
                ggplot2::theme_classic() +
                ggplot2::theme(
                  axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
                )

              # Case 2: Condition provided, plot grouped means with stats
            } else {
              if (!condition %in% names(colData(object))) { stop("'", condition, "' not found in colData slot.") }

              condition_df <- colData(object) %>%
                as.data.frame() %>%
                dplyr::mutate(SampleID = as.character(rownames(.)))

              plotting_df <- dplyr::left_join(plotting_df, condition_df, by = "SampleID")

              # --- NEW LOGIC: Filter by selected_groups if provided ---
              if (!is.null(selected_groups)) {

                # 1. Check which groups exist in the data
                available_groups <- unique(plotting_df[[condition]])

                # 2. Check for groups requested that don't exist
                missing_groups <- setdiff(selected_groups, available_groups)

                if (length(missing_groups) > 0) {
                  warning(paste0("The following items in 'selected_groups' were not found in the '",
                                 condition, "' column: ", paste(missing_groups, collapse = ", ")))
                }

                # 3. Filter the dataframe
                plotting_df <- plotting_df %>%
                  dplyr::filter(!!rlang::sym(condition) %in% selected_groups)

                if(nrow(plotting_df) == 0) {
                  stop("No data remaining after filtering by 'selected_groups'.")
                }
              }
              # --------------------------------------------------------

              groups <- unique(na.omit(plotting_df[[condition]]))

              # Main plot call
              g <- ggplot2::ggplot(plotting_df,
                                   ggplot2::aes(x = !!rlang::sym(condition), y = Intensity, fill = !!rlang::sym(condition))
              ) +
                # Draw the bars FIRST
                ggplot2::stat_summary(fun = mean, geom = "bar", color = "black") +
                # Draw ALL error bars SECOND, so they are on top.
                ggplot2::stat_summary(
                  fun.data = ggplot2::mean_sdl, fun.args = list(mult = 1),
                  geom = "errorbar", width = 0.2, color = "black"
                ) +
                ggplot2::labs(
                  title = paste("Mean Intensity for Protein:", prot),
                  subtitle = paste("Grouped by:", condition),
                  x = condition,
                  y = "Mean Intensity (+/- SD)"
                ) +
                ggplot2::theme_classic() +
                ggplot2::theme(legend.position = "none")

              # Pre-calculate stats and filter for ONLY significant pairs
              # We check if there are at least 2 groups REMAINING after filtering
              if (length(groups) >= 2) {
                stat_test <- tryCatch({

                  # --- This is the code that might fail ---
                  ggpubr::compare_means(as.formula(paste("Intensity ~", condition)), data = plotting_df) %>%
                    dplyr::filter(p.signif != "ns")
                  # ----------------------------------------

                }, error = function(e) {

                  # If an error occurs, stop the execution and provide a clear message.
                  stop(paste("Statistical comparison failed. This can happen if a group has too few non-missing observations.\nOriginal error:", e$message))

                })

                # Only add the layer if there are significant results to show
                if (nrow(stat_test) > 0) {
                  # Convert the significant pairs data.frame into the list format ggpubr requires
                  significant_comparisons <- lapply(1:nrow(stat_test), function(i) {
                    c(stat_test$group1[i], stat_test$group2[i])
                  })

                  g <- g + ggpubr::stat_compare_means(
                    comparisons = significant_comparisons,
                    label = "p.signif"
                  )
                }
              }
            }

            return(g)
          })
