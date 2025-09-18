#' Add a Processing Step to an Object's History
#'
#' Appends a formatted log entry to the `processing` slot of an S4 object.
#' This function is intended to be used internally by preprocessing functions
#' (e.g., `log2_transform`, `impute_data`) to record their actions.
#'
#' @param object An S4 object that contains a `processing` slot (list).
#' @param log_entry A named list containing metadata about the processing step.
#'   This list should ideally contain elements like `name`, `timestamp`, and
#'   `parameters`.
#'
#' @return The modified S4 object with the `log_entry` appended to its
#'   `processing` history.
#'
#' @keywords internal
add_processing_step <- function(object, log_entry) {
  object@processing <- append(object@processing, list(log_entry))
  return(object)
}

### Filtering #################################################################

#' @title Filter Proteins by Percentage of Valid Values
#'
#' @description
#' A generic function to filter proteins (rows) from a proteomics data object
#' based on the percentage of samples where they have valid observations.
#'
#' @param object A data object containing protein expression data.
#' @param percent The minimum percentage (from 0 to 100) of samples in which a
#'   protein must have a valid (non-NA) value to be retained.
#'
#' @return A modified object of the same class with proteins filtered.
#'
#' @export
setGeneric("filter_proteins_by_percent",
           def = function(object, percent) {
             standardGeneric("filter_proteins_by_percent")
           }
)

setMethod("filter_proteins_by_percent",
          signature(object = "ProtData", percent = "numeric"),
          function(object, percent) {

            # --- 1. Input Validation ---
            if (length(percent) != 1 || !is.numeric(percent) || percent < 0 || percent > 100) {
              stop("'percent' must be a single numeric value between 0 and 100.", call. = FALSE)
            }

            initial_protein_count <- nrow(object@data)
            if (initial_protein_count == 0) {
              warning("The data slot has no proteins to filter.", call. = FALSE)
              return(object)
            }

            # --- 2. Core Filtering Logic ---
            num_samples <- ncol(object@data)

            # Calculate the absolute number of samples required, taking the ceiling
            min_samples_required <- ceiling((percent / 100) * num_samples)

            message(paste0("Filtering to keep proteins with valid values in at least ", percent,
                           "% of samples (", min_samples_required, " out of ", num_samples, " samples)."))

            # Calculate the number of non-NA values for each row
            valid_counts <- rowSums(!is.na(object@data))

            # Create a logical vector indicating which proteins to keep
            proteins_to_keep <- valid_counts >= min_samples_required
            num_removed <- initial_protein_count - sum(proteins_to_keep)

            if (num_removed == 0) {
              message(paste0("All ", initial_protein_count, " proteins met the criterion. No proteins were removed."))
            } else {
              message(paste0(num_removed, " out of ", initial_protein_count, " proteins were removed."))
            }

            # --- 3. Subset the data and metadata slots ---
            object@data <- object@data[proteins_to_keep, , drop = FALSE]
            object@prot_meta <- object@prot_meta[proteins_to_keep, , drop = FALSE]

            # --- 4. Log the Operation ---
            log_entry <- list(
              name = "filter_proteins_by_percent",
              parameters = list(percent_threshold = percent,
                                calculated_min_samples = min_samples_required),
              details = paste(num_removed, "proteins removed.")
            )

            object <- add_processing_step(object, log_entry)

            # --- 5. Return the modified object ---
            return(object)
          }
)

#' Removes duplicate analytes
#'
#' Given a prot_meta column name, it will remove rows so that each value in the column is unique.
#' For each row the analyte with 1) the lowest missing values and 2) the greatest median intensity will be kept.
#'
#' @param PD A ProtData object.
#' @param col a column name of the prot_meta slot
#'
#' @return
#' @export
#'
setGeneric("filter_unique_proteins", function(object, col = NULL) standardGeneric("filter_unique_proteins"))
setMethod("filter_unique_proteins", "ProtData", function(object, col = NULL) {

  # --- Input Validation ---
  if (length(col) != 1) {
    stop("'col' must be a single character string.")
  }
  if (ncol(object@prot_meta)==0) {
    stop("The prot_meta slot is empty, there is no way to find duplicates")
  }
  if (!(col %in% names(object@prot_meta))) {
    stop(paste0("'", col, "' is not a valid column in the 'prot_meta' slot."))
  }

  if (is.null(col)) {
    col <- colnames(object@prot_meta)[1]
  }

  duplicates = as.numeric(nrow(object@prot_meta) - length(unique(object@prot_meta[[col]])))

  intensity_cols <- 1:ncol(object@data)

  # Combine data and metadata directly
  dat <- cbind(object@data, object@prot_meta) %>%
    as.data.frame() %>% # Ensure it's a dataframe for dplyr
    dplyr::mutate(missing_value = rowSums(is.na(dplyr::select(., all_of(intensity_cols))))) %>%
    dplyr::mutate(median = matrixStats::rowMedians(as.matrix(dplyr::select(., all_of(intensity_cols))), na.rm = TRUE)) %>%
    dplyr::group_by(!!as.name(col)) %>%
    dplyr::filter(missing_value == min(missing_value)) %>%
    dplyr::filter(median == max(median)) %>%
    dplyr::ungroup()

  # Separate back into data and metadata
  object@data <- dplyr::select(dat, all_of(intensity_cols))
  metadata_start_col <- length(intensity_cols) + 1
  metadata_end_col <- ncol(dat) - 2 # Subtract the two new columns (missing_value, median)
  object@prot_meta <- dplyr::select(dat, metadata_start_col:metadata_end_col)

  log_entry <- list(
    name = "filter_unique_proteins",
    parameters = list(condition = col),
    details = paste(duplicates, "proteins removed.")
  )

  object <- add_processing_step(object, log_entry)

  # --- 5. Return the modified object ---
  return(object)
})

#' Filter Outlier Samples Based on Protein Counts
#'
#' @description
#' This function identifies and removes entire samples that are considered outliers
#' based on the total number of non-missing proteins identified within them.
#'
#' The method calculates the mean and standard deviation of protein counts across
#' all samples. A sample is flagged as an outlier if its total protein count
#' falls outside the range defined by a specified number of standard deviations
#' (`sds`) from the mean.
#'
#' @param object A `ProtData` object.
#' @param sds The number of standard deviations from the mean protein count to use
#'   as the outlier threshold. Defaults to 3.
#'
#' @return A `ProtData` object with the outlier samples removed from the `data`
#'   and `condition` slots.
#'
#' @export
#'
#' @examples
#' # Create data where one sample has a very low number of identified proteins
#' raw_data <- data.frame(
#'   Gene = c("GENEA", "GENEB", "GENEC", "GENED"),
#'   SampleA = c(10, 11, 12, 13),
#'   SampleB = c(12, 13, 14, 15),
#'   SampleC = c(11, 12, 13, 14),
#'   SampleD_outlier = c(NA, NA, NA, 5) # Only one protein found
#' )
#'
#' pd_obj <- create_protdata(dat = raw_data)
#' cat("Dimensions before removing outliers:", dim(pd_obj@data), "\n")
#'
#' # With sds=1, the outlier sample should be identified and removed.
#' filtered_obj <- remove_outliers(pd_obj, sds = 1)
#' cat("Dimensions after removing outliers:", dim(filtered_obj@data), "\n")
#' print(colnames(filtered_obj@data))
#'
setGeneric("filter_outlier_samples", function(object, sds = 3) standardGeneric("filter_outlier_samples"))

setMethod("filter_outlier_samples",
          "ProtData",
          function(object, sds = 3){
            N_values <- colSums(!is.na(getData(object)))
            pgcounts <- data.frame(Sample = names(N_values), N = N_values)

            dat <- object@data
            cond <- object@condition

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
              dat <- dat[,!(colnames(dat) %in% outliers),drop = FALSE]
              cond <- cond[!rownames(cond) %in% outliers, ,drop = FALSE]

              object@data <- dat
              object@condition <- cond
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



#' get_overlap method for protdata class
#'
#' Filters the protdata object to retain only proteins that are present (non-NA)
#' in at least one sample within each unique group of the specified condition.
#'
#' @param object A protdata object.
#' @param condition_name A character string specifying the column name in the
#' condition' slot to group samples by.
#'
#' @return A new, modified protdata object containing only the overlapping proteins.
setGeneric("get_overlap",
           def = function(object, condition_name) {
             standardGeneric("get_overlap")
           }
)

setMethod("get_overlap",
          signature(object = "ProtData", condition_name = "character"),
          function(object, condition_name) {

            # --- Input Validation ---
            if (length(condition_name) != 1) {
              stop("'condition_name' must be a single character string.")
            }
            if (!(condition_name %in% names(object@condition))) {
              stop(paste0("'", condition_name, "' is not a valid column in the 'condition' slot."))
            }

            cat("Starting overlap analysis for condition: '", condition_name, "'\n", sep = "")

            # --- Extract data from slots ---
            cond_df <- object@condition
            data_df <- object@data
            original_prots <- nrow(data_df)

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

              # Get the names of the proteins that are present
              present_protein_names <- rownames(data_subset)[is_present]

              present_proteins_per_group[[as.character(group)]] <- present_protein_names

              cat("  - Group '", group, "': Found ", length(present_protein_names), " present proteins.\n", sep = "")
            }

            # --- Find the intersection of all sets of proteins ---
            if (length(present_proteins_per_group) == 0) {
              stop("No groups found for the specified condition.")
            }

            # Use Reduce with intersect to find common proteins across all lists
            overlapping_proteins <- Reduce(intersect, present_proteins_per_group)

            cat("Found", length(overlapping_proteins), "proteins present across all groups.\n")

            if (length(overlapping_proteins) == 0) {
              warning("No overlapping proteins found. Returning an empty object.")
            }

            # --- Modify the object slots based on the overlap ---
            # Note: In R, we return a modified copy rather than modifying in-place.
            # The user will typically assign the result back to the original variable.

            object@prot_meta <- object@prot_meta[overlapping_proteins, , drop = FALSE]
            object@data <- object@data[overlapping_proteins, , drop = FALSE]

            log_entry <- list(
              name = "get_overlap",
              parameters = list(condition = condition_name),
              details = paste(original_prots - length(overlapping_proteins), "proteins removed.")
            )

            object <- add_processing_step(object, log_entry)

            # --- 5. Return the modified object ---
            return(object)
          }
)
### Normalization #################################################################
#' Z-Score Normalization for Proteins Across Samples
#'
#' @description
#' This method performs a Z-score transformation on a protein-wise (row-wise)
#' basis. For each protein, it calculates the mean and standard deviation of its
#' abundance across all samples and then scales the values accordingly.
#'
#' This is a standard transformation for visualizing expression patterns in a
#' heatmap, as it highlights the relative change of each protein across samples,
#' independent of its absolute abundance.
#'
#' @param object A `ProtData` object.
#'
#' @return A `ProtData` object where the abundance values in the `@data` slot
#'   have been replaced by their row-wise Z-scores.
#'
#' @export
#' @rdname scale-ProtData
#' @aliases scale,ProtData-method
#'
#' @seealso [base::scale()]
#'
#' @examples
#' # Create sample data with proteins as rows
#' df <- data.frame(
#'   SampleA = c(100, 250, 50),
#'   SampleB = c(120, 200, 100),
#'   SampleC = c(110, 225, 75),
#'   row.names = c("Protein1", "Protein2", "Protein3")
#' )
#'
#' conditions <- data.frame(
#'   row.names = colnames(df),
#'   group = c("Control", "Treatment", "Control")
#' )
#'
#' prot_obj <- new("ProtData",
#'                 data = df,
#'                 condition = conditions,
#'                 method = "MS")
#'
#' # Check the means of each protein (row) before scaling
#' rowMeans(prot_obj@data)
#'
#' # Apply the scaling method
#' scaled_prot_obj <- scale(prot_obj)
#'
#' # The new data has row means near zero and row standard deviations of one
#' print(scaled_prot_obj@data)
#' cat("Row means after scaling:\n")
#' print(rowMeans(scaled_prot_obj@data))
#' cat("\nRow standard deviations after scaling:\n")
#' print(apply(scaled_prot_obj@data, 1, sd))
#'
setGeneric("z_score", function(object) standardGeneric("z_score"))

setMethod("z_score", "ProtData", function(object) {
  # Transpose, scale columns (which are now proteins), and transpose back
  object@data <- t(base::scale(t(object@data))) %>%
    as.data.frame()
  log_entry <- list(
    name = "z_score"
  )

  object <- add_processing_step(object, log_entry)

  # --- 5. Return the modified object ---
  return(object)
})

# normalization methods

#' Median Normalization of Proteomics Data
#'
#' @description
#' Performs median normalization on the quantitative data within a `ProtData`
#' object. This method corrects for systematic, sample-specific biases (e.g.,
#' differences in sample loading or instrument sensitivity) to make the samples
#' more comparable.
#'
#' The function operates by first calculating the median intensity for each
#' sample (column). It then calculates the median of these medians (the "global
#' median"). Finally, it scales the intensities in each sample by a
#' multiplicative factor so that every sample has the same median abundance,
#' equal to the global median.
#'
#' @details
#' This type of multiplicative adjustment is typically appropriate for raw,
#' non-log-transformed intensity data.
#'
#' @param object A `ProtData` object containing the abundance data to be
#'   normalized.
#'
#' @return A `ProtData` object with the abundance data in the `@data` slot
#'   normalized by the median-centering method.
#'
#' @export
#' @rdname median_normalize-ProtData
#' @aliases median_normalize,ProtData-method
#'
#'
setGeneric("median_normalize", function(object) standardGeneric("median_normalize"))
setMethod("median_normalize", "ProtData", function(object) {
  medians <- apply(object@data, 2, median, na.rm = TRUE) # per-sample medians
  global_median <- median(medians, na.rm = TRUE)
  object@data <- sweep(object@data, 2, medians, function(x, m) x * (global_median / m))

  log_entry <- list(
    name = "median_normalize",
    details = paste("global median: ", global_median)
  )

  object <- add_processing_step(object, log_entry)

  # --- 5. Return the modified object ---
  return(object)
})

#' Mean Normalization of Proteomics Data
#'
#' @description
#' Corrects for sample-specific biases by scaling each sample (column) to have
#' the same mean abundance.
#'
#' This method adjusts each sample's intensities by a multiplicative factor,
#' aligning the mean of every sample to a global mean calculated from the
#' entire dataset.
#'
#' @details
#' This multiplicative adjustment is typically used for raw, non-log-transformed
#' intensity data.
#'
#' @param object A `ProtData` object containing abundance data.
#'
#' @return A `ProtData` object with its data normalized by the mean-centering
#'   method.
#'
#' @export
#' @rdname mean_normalize-ProtData
#' @aliases mean_normalize,ProtData-method
#'
setGeneric("mean_normalize", function(object) standardGeneric("mean_normalize"))
setMethod("mean_normalize", "ProtData", function(object) {
  means <- apply(object@data, 2, mean, na.rm = TRUE)  # per-sample means
  global_mean <- mean(means, na.rm = TRUE)            # global mean
  object@data <- sweep(object@data, 2, means, function(x, m) x * (global_mean / m))
  log_entry <- list(
    name = "mean_normalize",
    details = paste("global mean: ", global_mean)
  )

  object <- add_processing_step(object, log_entry)

  # --- 5. Return the modified object ---
  return(object)
})

### Transformation #################################################################
#' Performs a log2 transform of protein intensity values
#'
#' @param PD A ProtData object.
#'
#' @return A ProtData object.
#' @export
#'
setGeneric("log2_transform", function(object) standardGeneric("log2_transform"))
setMethod("log2_transform", "ProtData", function(object) {
  object@data <- object@data %>%
    dplyr::mutate(across(where(is.numeric), ~ ifelse(!is.na(.) & !is.nan(.), log2(.+1), .)))
  log_entry <- list(
    name = "log2_transform"
  )

  object <- add_processing_step(object, log_entry)

  # --- 5. Return the modified object ---
  return(object)
})

#' Performs a log10 transform of protein intensity values
#'
#' @param PD A ProtData object.
#'
#' @return A ProtData object.
#' @export
#'
setGeneric("log_transform", function(object) standardGeneric("log_transform"))
setMethod("log_transform", "ProtData", function(object) {
  object@data <- object@data %>%
    dplyr::mutate(across(where(is.numeric), ~ ifelse(!is.na(.) & !is.nan(.), log(.+1), .)))
  log_entry <- list(
    name = "log10_transform"
  )

  object <- add_processing_step(object, log_entry)

  # --- 5. Return the modified object ---
  return(object)
})

### Imputation #################################################################
# imputation methods

#' Impute Missing Values with a Constant
#'
#' @description
#' Replaces all missing values (`NA` and `NaN`) in the numeric columns of the
#' data with a single, user-specified constant.
#'
#' @param object A `ProtData` object containing data with missing values.
#' @param value The numeric constant to use for replacing `NA` and `NaN` values.
#'
#' @return A `ProtData` object with missing values imputed.
#'
#' @export
#' @rdname impute-ProtData
#' @aliases impute,ProtData-method
#'
#' @examples
#' # Create a sample data frame with metadata and a missing numeric value
#' raw_data <- data.frame(
#'   Protein.ID = c("P02768", "P01023", "P60709"),
#'   Gene.Name = c("ALB", "A2M", "ACTB"),
#'   Sample_A = c(1.2e6, 2.3e6, NA),
#'   Sample_B = c(1.4e6, 2.6e6, 4.8e6)
#' )
#'
#' # Use the constructor to create a ProtData object.
#' # The constructor will automatically separate metadata from numeric data.
#' pd_obj <- create_protdata(dat = raw_data)
#'
#' # Impute the NA with 0
#' imputed_obj <- impute(pd_obj, value = 0)
#'
#' # View the imputed data slot
#' print(imputed_obj@data)
#'
setGeneric("impute", function(object, value) standardGeneric("impute"))
setMethod("impute", "ProtData", function(object, value) {
  total_missing <- sum(is.na(object@data))
  object@data <- object@data %>%
    dplyr::mutate(across(where(is.numeric), ~ ifelse(is.na(.) | is.nan(.), value, .)))
  log_entry <- list(
    name = "impute",
    parameters = list(fixed_value = value),
    details = paste(total_missing, "values imputed.")
  )

  object <- add_processing_step(object, log_entry)

  # --- 5. Return the modified object ---
  return(object)
})


#' Impute Missing Values with the Row Minimum
#'
#' @description
#' Performs row-wise minimum imputation. For each protein (row), it replaces
#' missing values (`NA`, `NaN`) with the minimum observed value found in that
#' same row.
#'
#' The imputation value can be scaled by a multiplicative factor `alpha`.
#'
#' @param object A `ProtData` object containing data with missing values.
#' @param alpha A numeric scaling factor to multiply the row minimum by before
#'   imputation. Defaults to 1 (no scaling).
#'
#' @return A `ProtData` object with missing values imputed on a per-protein basis.
#'
#' @export
#' @rdname impute_min-ProtData
#' @aliases impute_min,ProtData-method
#'
#' @examples
#' # Create data with different minimums and NAs in each row
#' raw_data <- data.frame(
#'   Gene = c("GENEA", "GENEB"),
#'   SampleA = c(100, 500),
#'   SampleB = c(200, 600),
#'   SampleC = c(NA, NA)
#' )
#'
#' pd_obj <- create_protdata(dat = raw_data)
#' cat("Original Data:\n")
#' print(pd_obj@data)
#'
#' # Impute using the row minimum (alpha = 1)
#' # Row 1's NA becomes 100; Row 2's NA becomes 500.
#' imputed_obj <- impute_min(pd_obj)
#' cat("\nImputed with alpha = 1:\n")
#' print(imputed_obj@data)
#'
#' # Impute using 90% of the row minimum
#' imputed_scaled <- impute_min(pd_obj, alpha = 0.9)
#' cat("\nImputed with alpha = 0.9:\n")
#' print(imputed_scaled@data)
#'
setGeneric("impute_min", function(object, alpha=1) standardGeneric("impute_min"))
setMethod("impute_min", "ProtData", function(object, alpha=1) {
  total_missing <- sum(is.na(object@data))
  object@data <- t(apply(object@data, 1, function(x) {
    min_val <- min(x[!is.na(x) & !is.nan(x)], na.rm = TRUE) * alpha  # find the minimum value excluding NA and NaN
    x[is.na(x) | is.nan(x)] <- min_val  # replace NA and NaN with the minimum value
    return(x)
  })) %>% as.data.frame()
  log_entry <- list(
    name = "impute_min",
    parameters = list(alpha = alpha),
    details = paste(total_missing, "values imputed.")
  )

  object <- add_processing_step(object, log_entry)

  # --- 5. Return the modified object ---
  return(object)
})


#' Impute from a Down-Shifted Normal Distribution
#'
#' @description
#' Performs row-wise imputation by drawing random values from a normal
#' distribution that is shifted to the left and narrower than the distribution of
#' observed values.
#'
#' This method assumes that missing values are primarily from proteins with low
#' abundance (i.e., below the detection limit). The default `shift` and `scale`
#' values are based on those used in the Perseus analysis platform.
#'
#' @details
#' This imputation method should be applied to log-transformed data, as the
#' underlying assumption of a normal distribution is more appropriate in log space.
#' An exception is made for rows with only one or zero observed values, where missing
#' values are imputed with 0.
#'
#' @param object A `ProtData` object containing data with missing values.
#' @param shift A numeric value specifying how many standard deviations to shift
#'   the mean of the distribution for imputed values. Default is 1.8.
#' @param scale A numeric value to scale the standard deviation of the
#'   distribution for imputed values. Default is 0.3.
#'
#' @return A `ProtData` object with missing values imputed from a simulated
#'   low-abundance distribution.
#'
#' @export
#' @rdname impute_left_dist-ProtData
#' @aliases impute_left_dist,ProtData-method
#'
#' @examples
#' # Create data with NAs, typically representing log-transformed values
#' raw_data <- data.frame(
#'   Gene = c("GENEA", "GENEB"),
#'   SampleA = c(25.1, 28.5),
#'   SampleB = c(25.5, 28.9),
#'   SampleC = c(NA, NA)
#' )
#'
#' pd_obj <- create_protdata(dat = raw_data)
#'
#' # For reproducibility of the random imputation
#' set.seed(123)
#'
#' imputed_obj <- impute_left_dist(pd_obj)
#' cat("Data after imputation:\n")
#' print(imputed_obj@data)
#'
setGeneric("impute_left_dist", function(object, shift = 1.8, scale = 0.3) standardGeneric("impute_left_dist"))
setMethod("impute_left_dist", "ProtData", function(object, shift = 1.8, scale = 0.3) {
  total_missing <- sum(is.na(object@data))
  imputed_data <- t(apply(object@data, 1, function(x) {
    mu <- mean(x, na.rm = TRUE)
    sigma <- sd(x, na.rm = TRUE)
    n_missing <- sum(is.na(x) | is.nan(x))

    if (n_missing == 0) {
      return(x)
    }
    else if (is.na(mu) || is.na(sigma)){
      x[is.na(x) | is.nan(x)] <- 0
      return(x)
    }

    # Draw random values, ensure they are non-negative
    imputed_vals <- rnorm(n_missing, mean = mu - shift * sigma, sd = sigma * scale)
    imputed_vals <- pmax(imputed_vals, 0)  # clip to zero

    x[is.na(x) | is.nan(x)] <- imputed_vals
    return(x)
  }))

  imputed_data <- as.data.frame(imputed_data)
  rownames(imputed_data) <- rownames(object@data)
  colnames(imputed_data) <- colnames(object@data)

  object@data <- imputed_data
  log_entry <- list(
    name = "impute_left_dist",
    parameters = list(shift = shift, scale = scale),
    details = paste(total_missing, "values imputed.")
  )

  object <- add_processing_step(object, log_entry)

  # --- 5. Return the modified object ---
  return(object)
})

### Correction #################################################################
#' Perform Batch Correction on a ProtData Object
#'
#' @description
#' This function adjusts for batch effects in the @data slot of a ProtData
#' object using the removeBatchEffect function from the limma package. It can
#' optionally preserve specified biological variation while removing the unwanted
#' technical variation associated with the batch variable.
#'
#' @param object An S4 object of class ProtData. The object must contain a
#'   @data slot with numeric data and a @condition slot with sample metadata.
#' @param batch_variable A single character string specifying the column name in
#'   the object's @condition slot that identifies the batch for each sample.
#' @param bio_variables An optional character vector of column names in the
#'   @condition slot that represent biological variables of interest. The
#'   variation from these variables will be preserved during the correction.
#'   If NULL (the default), a simpler correction is performed assuming a common
#'   mean across all samples.
#'
#' @return The input ProtData object with its @data slot updated with the
#'   batch-corrected values.
#'
#' @export
#'
#' @importFrom limma removeBatchEffect
#'
setGeneric("batch_correct",
           def = function(object, batch_variable, bio_variables = NULL) {
             standardGeneric("batch_correct")
           }
)

# --- Method 1: WHEN bio_variables ARE PROVIDED (No change) ---
setMethod("batch_correct",
          # The "ANY" signature allows this method to handle all cases.
          signature(object = "ProtData", batch_variable = "character", bio_variables = "ANY"),
          function(object, batch_variable, bio_variables) {
            # --- Input Validation for required arguments ---
            if (!requireNamespace("limma", quietly = TRUE)) stop("Please install the 'limma' package from Bioconductor.")
            if (length(batch_variable) != 1) stop("'batch_variable' must be a single column name.")
            if (!batch_variable %in% names(object@condition)) stop("'", batch_variable, "' not found in the condition slot.")
            if (!is.numeric(as.matrix(object@data))) stop("The data slot must be entirely numeric.")

            data_matrix <- as.matrix(object@data)
            batch_vector <- object@condition[[batch_variable]]

            # --- Use if/else to handle the optional 'bio_variables' ---

            # CASE 1: bio_variables WERE provided.
            if (!is.null(bio_variables)) {

              cat("Saatarting batch correction...\n")
              cat("  - Preserving biological variables:", paste(bio_variables, collapse=", "), "\n")

              # Validate the provided bio_variables
              if (!all(bio_variables %in% names(object@condition))) {
                stop("One or more 'bio_variables' not found in the condition slot.")
              }

              # Build the biological design matrix
              formula_str <- paste("~", paste(bio_variables, collapse = " + "))
              design_matrix <- model.matrix(as.formula(formula_str), data = object@condition)

              # Call limma WITH the design matrix
              corrected_data <- limma::removeBatchEffect(
                x = data_matrix,
                batch = batch_vector,
                design = design_matrix
              )

              # CASE 2: bio_variables were NOT provided (it is NULL).
            } else {

              cat("Starting batch correction...\n")
              cat("  - No biological variables provided. Relying on limma's default (intercept-only) design.\n")

              # Call limma WITHOUT the design matrix, relying on its default
              corrected_data <- limma::removeBatchEffect(
                x = data_matrix,
                batch = batch_vector
              )
            }

            # --- Update and Return the Object (same for both cases) ---
            object@data <- as.data.frame(corrected_data)
            cat("Batch correction complete. The data slot has been updated.\n")
            log_entry <- list(
              name = "batch_corrected",
              parameters = list(batch_variable = batch_variable, bio_variables =  bio_variables)
            )

            object <- add_processing_step(object, log_entry)

            # --- 5. Return the modified object ---
            return(object)
          }
)
