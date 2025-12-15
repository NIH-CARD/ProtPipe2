#' Add a Processing Step to an Object's History
#'
#' Appends a formatted log entry to the `processing` slot of an S4 object.
#' This function is intended to be used internally by preprocessing functions
#' (e.g., `log2_transform`, `impute_data`) to record their actions.
#'
#' @param object A \code{SummarizedExperiment} object with a processing_log entry in the metadata slot.
#' @param log_entry A named list containing metadata about the processing step.
#'   This list should ideally contain elements like `name`, `description`, and
#'   `parameters`.
#'
#' @return The modified object with the `log_entry` appended to its
#'   `processing_log` history.
#'
#' @keywords internal
add_processing_step <- function(object, log_entry) {
  metadata(object)$processing_log <- append(metadata(object)$processing_log, list(log_entry))
  return(object)
}

#' Generate Markdown Report of Preprocessing Steps
#'
#' @description
#' Takes a SummarizedExperiment object with a populated `processing_log`
#' in its metadata and writes a Markdown file summarizing the preprocessing
#' steps in order.
#'
#' @param object A SummarizedExperiment object.
#' @param output_file A character string giving the name of the output Markdown file.
#'
#' @return Invisibly returns the file path to the Markdown report.
#' @export
generate_preprocessing_report <- function(object, output_file = "preprocessing_report.md") {
  if (!"SummarizedExperiment" %in% class(object)) {
    stop("Input must be a SummarizedExperiment object.")
  }

  log <- metadata(object)$processing_log

  if (is.null(log) || length(log) == 0) {
    warning("No processing log found in metadata(object)$processing_log.")
    cat("# Preprocessing Report\n\nNo preprocessing steps recorded.\n", file = output_file)
    return(invisible(output_file))
  }

  # --- Build Markdown content ---
  lines <- c("# Preprocessing Report\n")
  lines <- c(lines, paste0("_Generated on ", Sys.Date(), "_\n\n"))

  for (i in seq_along(log)) {
    step <- log[[i]]

    # Step title
    step_title <- paste0("### Step ", i, ": ", step$name, "\n\n")

    # Parameters (formatted as list)
    if (!is.null(step$parameters) && length(step$parameters) > 0) {
      param_lines <- paste(
        names(step$parameters),
        ": ",
        unlist(step$parameters),
        collapse = "\n\n- "
      )
      param_lines <- paste0("- ", param_lines)
      param_lines <- c("**Parameters:\n\n", param_lines)
    } else {
      param_lines <- ""
    }

    if (!is.null(step$details) && length(step$details) > 0) {
      detail_lines <- paste(
        unlist(step$details),
        collapse = "\n\n- "
      )
      detail_lines <- paste0("- ", detail_lines)
      detail_lines <- c("**Details:\n\n", detail_lines)
    } else {
      detail_lines <- ""
    }

    # Details
    #details <- if (!is.null(step$details)) paste0( "**Details:\n\n", step$details) else ""

    lines <- c(lines,
               step_title,
               param_lines,
               detail_lines)
  }

  # --- Write Markdown file ---
  cat(paste(lines, collapse = "\n\n"), file = output_file)

  message("Preprocessing report written to: ", output_file)
  invisible(output_file)
}

### minimum intensity filtering #################################################################

#' @describeIn lod_filter Method for SummarizedExperiment objects
#' @export
setMethod("lod_filter", signature(se = "SummarizedExperiment"),
          function(se, lod_col = "Buffer") {

            # 1. Input Validation
            # (Note: The S4 signature handles class checking, but we keep your logic as requested)
            if (!is(se, "SummarizedExperiment")) {
              stop("Input 'se' must be a SummarizedExperiment object.")
            }

            lod_values <- NULL
            source_found <- "none"

            # 2. Search in rowData (metadata for rows/proteins)
            if (lod_col %in% colnames(rowData(se))) {
              lod_values <- rowData(se)[[lod_col]]
              source_found <- "rowData"
            }
            # 3. Search in Assay Columns (specific sample acting as LOD, e.g., Buffer)
            else if (lod_col %in% colnames(se)) {
              lod_values <- assay(se, 1)[, lod_col]
              source_found <- "assay"
            }
            else {
              stop(paste("Could not find", lod_col, "in rowData or assay columns of the SummarizedExperiment."))
            }

            # Check if lod_values is numeric
            if (!is.numeric(lod_values)) {
              stop(paste("The values in", lod_col, "are not numeric and cannot be used for filtering."))
            }

            # 4. Apply Filter to all assays
            # We iterate through all assays in the object (e.g., counts, log-intensities)
            assay_names <- assayNames(se)
            if (is.null(assay_names)) assay_names <- paste0("assay_", seq_along(assays(se)))

            mat <- assay(se)
            mask <- mat < lod_values
            num_low <- sum(mask, na.rm = TRUE)
            mask[is.na(mask)] <- FALSE
            mat[mask] <- NA
            assay(se) <- mat

            # Log the Operation ---
            log_entry <- list(
              name = "lod_filter",
              details = paste(num_low, "proteins removed.")
            )

            object <- add_processing_step(se, log_entry)
            return(object)
          }
)

#' @describeIn apply_min_intenisty Method for SummarizedExperiment objects
#' @export
setMethod("apply_min_intenisty", signature(object = "SummarizedExperiment"),
          function(object, lod) {

            mat <- assay(object)
            num_low <- sum(mat < lod, na.rm = TRUE)
            mat[mat < lod] <- NA
            assay(object) <- mat

            # Log the Operation ---
            log_entry <- list(
              name = "apply_min_intenisty",
              parameters = list(LOD = lod),
              details = paste(num_low, "proteins removed.")
            )

            object <- add_processing_step(object, log_entry)
            return(object)
          }
)

### outlier removal #################################################################

#' @describeIn filter_proteins_by_percent
#' @export
setMethod("filter_proteins_by_percent",
          signature(object = "SummarizedExperiment", percent = "numeric"),
          function(object, percent) {

            # --- 1. Input Validation ---
            if (length(percent) != 1 || !is.numeric(percent) || percent < 0 || percent > 100) {
              stop("'percent' must be a single numeric value between 0 and 100.", call. = FALSE)
            }

            initial_protein_count <- nrow(assay(object))
            if (initial_protein_count == 0) {
              warning("The data slot has no proteins to filter.", call. = FALSE)
              return(object)
            }

            # --- 2. Core Filtering Logic ---
            num_samples <- ncol(assay(object))

            # Calculate the absolute number of samples required, taking the ceiling
            min_samples_required <- ceiling((percent / 100) * num_samples)

            message(paste0("Filtering to keep proteins with valid values in at least ", percent,
                           "% of samples (", min_samples_required, " out of ", num_samples, " samples)."))

            # Calculate the number of non-NA values for each row
            valid_counts <- rowSums(!is.na(assay(object)))

            # Create a logical vector indicating which proteins to keep
            proteins_to_keep <- valid_counts >= min_samples_required
            num_removed <- initial_protein_count - sum(proteins_to_keep)

            if (num_removed == 0) {
              message(paste0("All ", initial_protein_count, " proteins met the criterion. No proteins were removed."))
            } else {
              message(paste0(num_removed, " out of ", initial_protein_count, " proteins were removed."))
            }

            # --- 3. Subset the data and metadata slots ---
            object <- object[proteins_to_keep, ]

            # --- 4. Log the Operation ---
            log_entry <- list(
              name = "filter_proteins_by_percent",
              parameters = list(percent_threshold = percent),
              details = paste(num_removed, "proteins removed.")
            )

            object <- add_processing_step(object, log_entry)

            # --- 5. Return the modified object ---
            return(object)
          }
)

#' @describeIn filter_unique_proteins
#' @export
setMethod("filter_unique_proteins",
          signature(object = "SummarizedExperiment", col = "character"),
          function(object, col = NULL) {

  # --- Input Validation ---
  if (length(col) != 1) {
    stop("'col' must be a single character string.")
  }
  if (ncol(rowData(object))==0) {
    stop("The rowData slot is empty, there is no way to find duplicates")
  }
  if (!(col %in% names(rowData(object)))) {
    stop(paste0("'", col, "' is not a valid column in the 'rowData' slot."))
  }

  if (is.null(col)) {
    col <- colnames(object@prot_meta)[1]
  }

  assay_data <- assay(object)
  missing_values <- rowSums(is.na(assay_data))
  median_intensities <- matrixStats::rowMedians(as.matrix(assay_data), na.rm = TRUE)

  # --- 2. Add Metrics to a Temporary rowData Copy ---
  # This avoids modifying the original object until the final step
  temp_rd <- rowData(object)
  temp_rd$missing_value <- missing_values
  temp_rd$median_intensity <- median_intensities
  temp_rd$original_index <- 1:nrow(temp_rd) # Keep track of original position

  # --- 3. Determine Indices of Rows to Keep ---
  # Use dplyr on the rowData to perform the filtering logic
  rows_to_keep_df <- as.data.frame(temp_rd) %>%
    dplyr::group_by(!!as.name(col)) %>%
    # First, filter for the minimum number of missing values within the group
    dplyr::filter(missing_value == min(missing_value)) %>%
    # Then, as a tie-breaker, filter for the highest median intensity
    dplyr::filter(median_intensity == max(median_intensity)) %>%
    # In case of a perfect tie, just take the first one
    dplyr::slice(1) %>%
    dplyr::ungroup()

  indices_to_keep <- rows_to_keep_df$original_index

  filtered_object <- object[indices_to_keep, ]

  duplicates <- nrow(object) - length(indices_to_keep)
  log_entry <- list(
    name = "filter_unique_proteins",
    parameters = list(condition = col),
    details = paste(duplicates, "proteins removed.")
  )

  filtered_object <- add_processing_step(filtered_object, log_entry)
  return(filtered_object)
})

#' @describeIn filter_outlier_samples
#' @export
setMethod("filter_outlier_samples",
          signature(object = "SummarizedExperiment", sds = "numeric"),
          function(object, sds = 3){
            N_values <- colSums(!is.na(assay(object)))
            pgcounts <- data.frame(Sample = names(N_values), N = N_values)

            dat <- assay(object)
            cond <- colData(object)

            stdev <- sd(pgcounts[,'N'])
            mean_count <- mean(pgcounts[,'N'])
            min_protein_groups <- floor(mean_count - (sds * stdev))
            max_protein_groups <- ceiling(mean_count + (sds * stdev))
            cat(paste0('INFO: Tolerating protein group counts in the range [', min_protein_groups,',',max_protein_groups,']\n'))
            low_count_samples <- as.character(pgcounts[pgcounts$N < min_protein_groups, 'Sample'])
            if(length(low_count_samples)>0){
              print("removing the following samples with low protein counts:")
              print(low_count_samples)
            }
            high_count_samples <- as.character(pgcounts[pgcounts$N > max_protein_groups, 'Sample'])
            if(length(high_count_samples)>0){
              print("removing the following samples with high protein counts:")
              print(high_count_samples)
            }

            outliers <- c(low_count_samples, high_count_samples)
            if (length(outliers) >0){
              object <- object[, !(colnames(dat) %in% outliers)]
            }
            log_entry <- list(
              name = "filter_outlier_samples",
              parameters = list(standard_deviations = sds),
              details = paste(length(outliers), "samples removed.")
            )

            object <- add_processing_step(object, log_entry)

            # --- 5. Return the modified object ---
            return(object)
          }
)

#' @describeIn filter_overlap
#' @export
setMethod("filter_overlap",
          signature(object = "SummarizedExperiment", condition_name = "character"),
          function(object, condition_name) {

            # --- Input Validation ---
            if (length(condition_name) != 1) {
              stop("'condition_name' must be a single character string.")
            }
            if (!(condition_name %in% names(colData(object)))) {
              stop(paste0("'", condition_name, "' is not a valid column in the 'condition' slot."))
            }

            cat("Starting overlap analysis for condition: '", condition_name, "'\n", sep = "")

            # --- Extract data from slots ---
            cond_df <- colData(object)
            data_df <- assay(object)
            original_prots <- nrow(object)

            # Get the unique groups from the specified condition column
            groups <- unique(cond_df[[condition_name]])

            # A list to store the names of present proteins for each group
            present_proteins_per_group <- list()

            # --- Logic to find present proteins for each group ---
            for (group in groups) {
              # Get the sample names belonging to the current group
              samples_in_group <- rownames(cond_df)[cond_df[[condition_name]] == group]

              # Subset the data to include only these samples
              # Use drop=FALSE to ensure it remains a data.frame even with one sample
              data_subset <- data_df[, samples_in_group, drop = FALSE]

              # A protein is "present" if it has no NA values across all samples in this group
              # This is a strict definition of presence.
              is_present <- rowSums(!is.na(data_subset)) > 0

              # # Get the names of the proteins that are present
              # present_protein_names <- rownames(data_subset)[is_present]

              present_proteins_per_group[[as.character(group)]] <- is_present#present_protein_names

              cat("  - Group '", group, "': Found ", sum(is_present), " present proteins.\n", sep = "")
            }

            # --- Find the intersection of all sets of proteins ---
            if (length(present_proteins_per_group) == 0) {
              stop("No groups found for the specified condition.")
            }

            # Use Reduce with intersect to find common proteins across all lists
            overlapping_proteins <- Reduce(`&`, present_proteins_per_group)

            cat("Found", sum(overlapping_proteins), "proteins present across all groups.\n")

            if (sum(overlapping_proteins) == 0) {
              warning("No overlapping proteins found. Returning an empty object.")
            }

            # --- Modify the object slots based on the overlap ---
            # Note: In R, we return a modified copy rather than modifying in-place.
            # The user will typically assign the result back to the original variable.
            object <- object[overlapping_proteins,]

            log_entry <- list(
              name = "filter_overlap",
              parameters = list(condition = condition_name),
              details = paste(original_prots - sum(overlapping_proteins), "proteins removed.")
            )

            object <- add_processing_step(object, log_entry)

            # --- 5. Return the modified object ---
            return(object)
          }
)
### Normalization #################################################################

#' @describeIn z_score
#' @export
setMethod("z_score",
          signature(object = "SummarizedExperiment"),
          function(object) {
  # Transpose, scale columns (which are now proteins), and transpose back
  assay(object) <- t(base::scale(t(assay(object))))
  log_entry <- list(
    name = "z_score"
  )

  object <- add_processing_step(object, log_entry)

  # --- 5. Return the modified object ---
  return(object)
})

# normalization methods


#' @describeIn median_normalize
#' @export
setMethod("median_normalize",
          signature(object = "SummarizedExperiment"),
          function(object) {
  medians <- apply(assay(object), 2, median, na.rm = TRUE) # per-sample medians
  global_median <- median(medians, na.rm = TRUE)
  assay(object) <- sweep(assay(object), 2, medians, function(x, m) x * (global_median / m))

  log_entry <- list(
    name = "median_normalize",
    details = paste("global median: ", global_median)
  )

  object <- add_processing_step(object, log_entry)

  # --- 5. Return the modified object ---
  return(object)
})


#' @describeIn "mean_normalize"
#' @export
setMethod("mean_normalize",
          signature(object = "SummarizedExperiment"),
          function(object) {
  means <- apply(assay(object), 2, mean, na.rm = TRUE)  # per-sample means
  global_mean <- mean(means, na.rm = TRUE)            # global mean
  assay(object) <- sweep(assay(object), 2, means, function(x, m) x * (global_mean / m))
  log_entry <- list(
    name = "mean_normalize",
    details = paste("global mean: ", global_mean)
  )

  object <- add_processing_step(object, log_entry)

  # --- 5. Return the modified object ---
  return(object)
})

### Transformation #################################################################

#' @describeIn log2_transform
#' @export
setMethod("log2_transform",
          signature(object = "SummarizedExperiment"),
          function(object) {
  assay(object) <- log2(assay(object) +1)

  log_entry <- list(
    name = "log2_transform"
  )

  object <- add_processing_step(object, log_entry)

  # --- 5. Return the modified object ---
  return(object)
})


### Imputation #################################################################


#' @describeIn impute
#' @export
setMethod(
  "impute",
  signature(object = "SummarizedExperiment", value = "numeric"),
  function(object, value) {
    stopifnot(length(value) == 1, is.numeric(value))
    m <- assay(object)
    m[m < 0] <- NA

    # NA *and* NaN are TRUE for is.na()
    total_missing <- sum(is.na(m))

    if (methods::is(m, "sparseMatrix")) {
      # Only stored entries are in @x; structural zeros stay zero
      m@x[is.na(m@x)] <- value
    } else {
      m[is.na(m)] <- value
    }

    assay(object) <- m

    log_entry <- list(
      name = "impute",
      parameters = list(fixed_value = value),
      details = paste(total_missing, "values imputed.")
    )
    object <- add_processing_step(object, log_entry)

    object
  }
)

#' @describeIn impute_min
#' @export
## Optional: use these if you have DelayedArray/DelayedMatrixStats installed
## library(DelayedArray)
## library(DelayedMatrixStats)
## library(matrixStats)

setMethod(
  "impute_min","SummarizedExperiment",
  function(object, alpha = 1) {
    stopifnot(length(alpha) == 1, is.numeric(alpha))

    m <- assay(object)
    m[m < 0] <- NA

    # total NA/NaN count (is.na() treats NaN as TRUE)
    total_missing <- sum(is.na(m))

    # Row minima excluding NAs
    row_mins <- if (methods::is(m, "DelayedMatrix")) {
      # Efficient for HDF5-backed/DelayedArray
      DelayedMatrixStats::rowMins(m, na.rm = TRUE)
    } else if (is.matrix(m)) {
      # Fast for in-memory matrices
      matrixStats::rowMins(m, na.rm = TRUE)
    } else {
      # Fallback
      apply(m, 1, function(x) min(x, na.rm = TRUE))
    }

    # Scale by alpha; rows that were all-NA have +Inf — don't impute those
    impute_vals <- alpha * row_mins
    bad_rows <- !is.finite(impute_vals)  # Inf (all-NA rows) or NA
    if (any(bad_rows)) impute_vals[bad_rows] <- NA_real_

    # Replace NA/NaN per row without changing dimnames
    miss_idx <- which(is.na(m), arr.ind = TRUE)
    if (nrow(miss_idx) > 0) {
      m[miss_idx] <- impute_vals[miss_idx[,"row"]]
    }

    assay(object) <- m

    log_entry <- list(
      name = "impute_min",
      parameters = list(alpha = alpha),
      details = paste(total_missing, "values imputed.")
    )
    object <- add_processing_step(object, log_entry)

    object
  }
)

#' @describeIn impute_left_dist
#' @export
setMethod(
  "impute_left_dist","SummarizedExperiment",
  function(object, shift = 1.8, scale = 0.3) {
    stopifnot(length(shift) == 1, length(scale) == 1,
              is.numeric(shift), is.numeric(scale))

    m <- assay(object)
    m[m < 0] <- NA
    total_missing <- sum(is.na(m))  # counts NA and NaN

    # Row-wise imputation without changing dimnames/class of 'm'
    for (i in seq_len(nrow(m))) {
      x <- m[i, ]
      miss <- is.na(x)
      n_missing <- sum(miss)
      if (n_missing == 0) next

      mu <- mean(x, na.rm = TRUE)
      sigma <- stats::sd(x, na.rm = TRUE)

      if (!is.finite(mu) || !is.finite(sigma)) {
        # All-NA (or degenerate) row: fill with 0
        imp <- rep(0, n_missing)
      } else {
        # Left-shifted normal; clip at 0
        imp <- pmax(stats::rnorm(n_missing, mean = mu - shift * sigma,
                                 sd = sigma * scale), 0)
      }
      m[i, miss] <- imp
    }

    assay(object) <- m

    log_entry <- list(
      name = "impute_left_dist",
      parameters = list(shift = shift, scale = scale),
      details = paste(total_missing, "values imputed.")
    )
    object <- add_processing_step(object, log_entry)

    object
  }
)

### Correction #################################################################

#' @describeIn batch_correct Method for SummarizedExperiment objects.
#' @export
setMethod("batch_correct",
          signature(object = "SummarizedExperiment", batch_variable = "character", bio_variables = "ANY"),
          function(object, batch_variable, bio_variables) {
            # --- 1. Input Validation ---
            if (!requireNamespace("limma", quietly = TRUE)) {
              stop("Please install the 'limma' package from Bioconductor.")
            }
            if (length(batch_variable) != 1) {
              stop("'batch_variable' must be a single column name.")
            }

            col_data <- colData(object)
            assay_data <- assay(object) # Use the first assay by default

            if (!batch_variable %in% colnames(col_data)) {
              stop("'", batch_variable, "' not found in colData.")
            }
            if (!is.numeric(assay_data)) {
              stop("The assay data must be numeric for batch correction.")
            }

            batch_vector <- col_data[[batch_variable]]

            # --- 2. Perform Batch Correction (using limma) ---

            # CASE 1: Biological variables WERE provided.
            if (!is.null(bio_variables)) {
              if (!is.character(bio_variables)) stop("'bio_variables' must be a character vector of column names.")
              if (!all(bio_variables %in% colnames(col_data))) {
                stop("One or more 'bio_variables' not found in colData.")
              }

              message("Starting batch correction...")
              message("  - Preserving biological variables: ", paste(bio_variables, collapse=", "))

              # Coerce to data.frame for model.matrix compatibility
              design_df <- as.data.frame(col_data)
              formula_str <- paste("~", paste(bio_variables, collapse = " + "))
              design_matrix <- model.matrix(as.formula(formula_str), data = design_df)

              corrected_data <- limma::removeBatchEffect(
                x = assay_data,
                batch = batch_vector,
                design = design_matrix
              )

              # CASE 2: Biological variables were NOT provided.
            } else {
              message("Starting batch correction...")
              message("  - No biological variables provided; protecting the intercept.")

              corrected_data <- limma::removeBatchEffect(
                x = assay_data,
                batch = batch_vector
              )
            }

            # --- 3. Update and Return the Object ---

            # Best practice: Store the corrected data in a NEW assay
            assay(object) <- corrected_data

            # Log the operation
            log_entry <- list(
              name = "batch_corrected",
              parameters = list(batch_variable = batch_variable, bio_variables = bio_variables)
            )

            # Assuming 'add_processing_step' is your helper that works with metadata()
            object <- add_processing_step(object, log_entry)

            return(object)
          }
)
