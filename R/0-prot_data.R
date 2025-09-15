#' @importFrom magrittr %>%
#' ProtData Class
#'
#' An S4 class that holds proteomics data and provides methods for processing.
#'
#' @export
setClass("ProtData",
         slots = list(
           data = "data.frame",      # The main proteomics data (proteins are rows, samples are columns)
           condition = "data.frame",  # Conditions of the samples. Rownames must match colnames of data
           prot_meta = "data.frame",  # information about the rows (genes, organism, etc...)
           method = "character"       # The method used (e.g., "MS", "Somascan", etc.)
         ),
         prototype = list(
           data = data.frame(),        # Default is an empty data frame
           condition = data.frame(),   # Default is an empty data frame
           prot_meta = data.frame(),   # Default is an empty data frame
           method = "Unknown"          # Default method is set to "Unknown"
         )
)


# Constructor for ProtData class
#' Create a ProtData Object
#'
#' This function creates an instance of the ProtData class.
#'
#' @param dat A data frame containing proteomics data (proteins are rows, samples are columns).
#' @intensity_cols vector of column indices corresponding to protein intensities. This will default to numeric columns.Optional.
#' @param condition A data frame containing conditions of the samples. Rownames should match colnames of data. Optional.
#' @param method A character string describing the method used for generating the data. Optional.
#'
#' @return An instance of the ProtData class.
#' @export
create_protdata <- function(dat, intensity_cols = NULL, condition = NULL, method = "Unknown") {



  # Check that data is a data frame
  if (!is.data.frame(dat)) {
    stop("The 'data' argument must be a data frame.")
  }

  # convert cols to numeric
  dat <- convert_numeric_cols(dat)

  #get the intensity cols
  if (is.null(intensity_cols)){
    intensity_cols <- detect_intensity_cols(dat)
  }


  #standardize the column names and order
  #dat <- standardize_format(data)
  colnames(dat) <- trim_names(colnames(dat))
  # col_order=c(colnames(dat)[1:2],sort(colnames(dat)[3:ncol(dat)]))
  # data.table::setcolorder(dat,col_order)

  # This should be a seperate function
  # dat <- dat %>%
  #   # Calculte missing values
  #   dplyr::mutate(missing_value = rowSums(is.na(dplyr::select(., -contains("Protein_Group|Genes")))))%>%
  #   # Calculate median values
  #   dplyr::mutate(median = matrixStats::rowMedians(as.matrix(dplyr::select(., -contains("Protein_Group|Genes|missing_value")) %>%
  #                                          dplyr::select_if(is.numeric)), na.rm = TRUE)) %>%
  #   #Identify unique Protein_Group with the least missing values and highest median intensity
  #   dplyr::group_by(Protein_Group) %>%
  #   dplyr::filter(missing_value == min(missing_value)) %>%
  #   dplyr::slice(which.max(median)) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::select(-c(missing_value, median))

  #Split up intensity and non-numeric protein metadata into seperate dataframes
  prot_meta <- as.data.frame(dat[, -intensity_cols])
  data <- as.data.frame(dat[, intensity_cols])
  data <- data %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.numeric(.)))

  ## CONDITION FILE ###################

  # Ensure that condition, if provided, has SampleID matching the colnames of data
  if (!is.null(condition)) {
    condition <- as.data.frame(condition)

    # make sample col the rownames if it exists
    if (!"SampleID" %in% colnames(condition)) {
      stop("Error: The condition file is missing the required 'SampleID' column.")
    }
    if (any(duplicated(condition$SampleID))) {
      stop("Error: The 'SampleID' column contains duplicate values. All sample IDs must be unique.")
    }
    rownames(condition) <- as.character(condition$SampleID)
    condition$SampleID <- NULL

    # trim the names
    rownames(condition) <- trim_names(rownames(condition))
    # drop excess rows from the condition file if they exist
    if (!all(rownames(condition) %in% colnames(data))) {
      print("Rownames of 'condition' do not match the colnames of 'data'.")

      # Drop rows of condition that are not in data
      matching_rows <- intersect(rownames(condition), colnames(data))
      condition <- condition[matching_rows, ]

      dropped_rows <- setdiff(rownames(condition), matching_rows)
      warning("Dropped rows:\n")
      warning(dropped_rows)
    }
    # add additional rows to condition if they exist in data
    if (ncol(data) > nrow(condition)) {
      # Find missing columns
      missing_cols <- setdiff(colnames(data), rownames(condition))

      # Create rows with NA for the missing columns and add them to 'condition'
      missing_rows <- data.frame(matrix(NA, nrow = length(missing_cols), ncol = ncol(condition)))
      colnames(missing_rows) <- names(condition)
      rownames(missing_rows) <- missing_cols

      # Add missing rows to 'condition'
      condition <- rbind(condition, missing_rows)
    }

    # Reorder 'condition' to match the order of 'data' columns
    condition <- condition[match(colnames(data), rownames(condition)), , drop = FALSE]
  }

  # if condition file is not provided, one is generated by removing _{replicate number} from each
  # sample name if present
  else {
    condition <- data.frame(
      base_condition = gsub("_\\d+$", "", colnames(data))  # Remove _ followed by digits at the end of the column names
    )
    rownames(condition) <- colnames(data)
  }


  # Create a new ProtData object
  new("ProtData",
      data = data,
      condition = condition,
      prot_meta = prot_meta,
      method = method)
}

#' Identify Numeric Columns by Index
#'
#' @description
#' A helper function that scans a data frame and returns the integer indices
#' of the columns that contain numeric data.
#'
#' This is primarily used internally by data loading functions to automatically
#' distinguish quantitative intensity columns from metadata columns.
#'
#' @param df A `data.frame` to be scanned.
#'
#' @return An integer vector containing the column indices of all numeric
#'   columns found in the input `data.frame`.
#'
#' @keywords internal
#'
#' @examples
#' # Create a sample data frame with mixed data types
#' sample_df <- data.frame(
#'   Protein = c("P02768", "P01023", "P10636"),
#'   Gene = c("ALB", "A2M", "ACTB"),
#'   Sample.1 = c(25.1, 22.4, 30.1),
#'   Sample.2 = c(26.2, 21.9, 31.5),
#'   IsContaminant = c(FALSE, FALSE, FALSE)
#' )
#'
#' # Use the function to find the numeric columns
#' # In this case, it should identify columns 3 and 4
#' ProtPipe:::detect_intensity_cols(sample_df)
#'
#' @export
detect_intensity_cols <- function(df) {
  which(sapply(df, is.numeric))
}



####### GETTERS and SETTERS ###############################################################

# Define the generic for 'getData'
setGeneric("getData", function(object) standardGeneric("getData"))

# Define the setter for 'getData'
setGeneric("setData", function(object, value) standardGeneric("setData"))

# Now define the method for 'getData' for the ProtData class
setMethod("getData", "ProtData", function(object) {
  return(object@data)
})

# Define the setter method for 'setData' for ProtData
setMethod("setData", "ProtData", function(object, value) {
  if (!inherits(value, "data.table")) {
    stop("'data' must be a data.table.")
  }
  object@data <- value
  return(object)
})

# Define the generic for 'data.long'
setGeneric("getDataLong", function(object) standardGeneric("getDataLong"))
setGeneric("setDataLong", function(object, value) standardGeneric("setDataLong"))

# Define the methods for 'data.long' for ProtData
setMethod("getDataLong", "ProtData", function(object) {
  return(object@data.long)
})

setMethod("setDataLong", "ProtData", function(object, value) {
  if (!inherits(value, "data.table")) {
    stop("'data.long' must be a data.table.")
  }
  object@data.long <- value
  return(object)
})

# Similarly, define getter and setter for 'condition'
setGeneric("getCondition", function(object) standardGeneric("getCondition"))
setGeneric("setCondition", function(object, value) standardGeneric("setCondition"))

setMethod("getCondition", "ProtData", function(object) {
  return(object@condition)
})

setMethod("setCondition", "ProtData", function(object, value) {
  if (!inherits(value, "data.frame")) {
    stop("'condition' must be a data.frame.")
  }
  object@condition <- value
  return(object)
})

# Similarly, define getter and setter for 'prot_meta'
setGeneric("getProtMeta", function(object) standardGeneric("getProtMeta"))
setGeneric("setProtMeta", function(object, value) standardGeneric("setProtMeta"))

setMethod("getProtMeta", "ProtData", function(object) {
  return(object@prot_meta)
})

setMethod("setProtMeta", "ProtData", function(object, value) {
  if (!inherits(value, "data.table")) {
    stop("'prot_meta' must be a data.table.")
  }
  object@prot_meta <- value
  return(object)
})

# Similarly, define getter and setter for 'method'
setGeneric("getProtMethod", function(object) standardGeneric("getProtMethod"))
setGeneric("setProtMethod", function(object, value) standardGeneric("setProtMethod"))

#WHY DOESNT THIS WORK????
setMethod("getProtMethod", "ProtData", function(object) {
  return(object@method)
})

setMethod("setProtMethod",
          signature = c(object = "ProtData", value = "character"),
          function(object, value) {
            if (length(value) != 1) {
              stop("'method' must be a single character string.")
            }
            object@method <- value
            return(object)
          })

####### Some Class Methods ###############################################################

#' Number of samples
#'
#' @param PD A ProtData object.
#'
#' @return number of samples
#' @export
#'
setGeneric("num_samples", function(object) standardGeneric("num_samples"))
setMethod("num_samples", "ProtData", function(object) {
  return(as.numeric(ncol(object@data)))
})

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
  return(object)
})

#' Z-Score Transformation for Proteins Across Samples
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
setGeneric("scale", function(object) standardGeneric("scale"))

setMethod("scale", "ProtData", function(object) {
  # Transpose, scale columns (which are now proteins), and transpose back
  object@data <- t(base::scale(t(object@data))) %>%
    as.data.frame()
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
  return(object)
})

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
  object@data <- object@data %>%
    dplyr::mutate(across(where(is.numeric), ~ ifelse(is.na(.) | is.nan(.), value, .)))
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
  object@data <- t(apply(object@data, 1, function(x) {
    min_val <- min(x[!is.na(x) & !is.nan(x)], na.rm = TRUE) * alpha  # find the minimum value excluding NA and NaN
    x[is.na(x) | is.nan(x)] <- min_val  # replace NA and NaN with the minimum value
    return(x)
  })) %>% as.data.frame()
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
  return(object)
})

# batch correct

# (Generic function definition remains the same)
#' Title
#'
#' @param object
#' @param batch_variable
#' @param bio_variables
#'
#' @return
#' @export
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
          ff<<-"kkas"
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
setGeneric("unique_data", function(object, col = NULL) standardGeneric("unique_data"))
setMethod("unique_data", "ProtData", function(object, col = NULL) {
  if (is.null(col)) {
    col <- colnames(object@prot_meta)[1]
  }
  intensity_cols <<- 1:ncol(object@data)

  dataa <<- object@data %>%
    as.data.frame() %>% # Ensure it's a dataframe for dplyr
    dplyr::mutate(missing_value = rowSums(is.na(dplyr::select(., all_of(intensity_cols))))) %>%
    dplyr::mutate(median = matrixStats::rowMedians(as.matrix(dplyr::select(., all_of(intensity_cols))), na.rm = TRUE))

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

  return(object)
})




######## HELPER FUNCTIONS ####################################################################

#total proteomics############
standardize_format <- function(DT.original) {
  # Accepts an input protein group intensity data.table, whether spectronaut or DIA-NN format,
  # and restructures into one consistent style for downstream processing
  DT <- DT.original
  if("Protein.Ids" %in% colnames(DT) | "Protein.Group" %in% colnames(DT)) {
    # print("DIAnn input")
    # DT[, 'Protein.Ids' := NULL]
    # DT[, 'Protein.Names' := NULL]
    # DT[, 'First.Protein.Description' := NULL]
    data.table::setnames(DT, 'Protein.Group', 'Protein_Group')
  }
  else if('EG.PrecursorId' %in% colnames(DT)) {
    print("Spectronaut input")
    data.table::setnames(DT, 'EG.PrecursorId', 'Peptide_Sequence')
    data.table::setnames(DT, 'PG.Genes', 'Genes')
    DT=as.data.frame(DT)
    # Use only Protein_Group and Genes
    dplyr::select=c('Peptide_Sequence','Genes',grep('raw',colnames(DT),value = T))
    DT=DT[, dplyr::select]
    #as number
    DT[,grep('raw',colnames(DT))]=as.data.frame(apply(DT[,grep('raw',colnames(DT))],2,as.numeric))
    DT=data.table(DT)
  }
  else if('PG.ProteinGroups' %in% colnames(DT)) {
    print("Spectronaut input")
    data.table::setnames(DT, 'PG.ProteinGroups', 'Protein_Group')
    data.table::setnames(DT, 'PG.Genes', 'Genes')
    DT=as.data.frame(DT)
    # Use only Protein_Group and Genes
    dplyr::select=c('Protein_Group','Genes',grep('PG.Quantity',colnames(DT),value = T))
    DT=DT[, dplyr::select]
    #as number
    DT[,grep('PG.Quantity',colnames(DT))]=as.data.frame(apply(DT[,grep('PG.Quantity',colnames(DT))],2,as.numeric))
    DT=data.table(DT)
  }
  else if('Peptide Sequence' %in% colnames(DT)) {
    print("FragPipe input")
    data.table::setnames(DT, 'Gene', 'Genes')
    colnames(DT)=gsub("\\s", "_",colnames(DT))
    # Use only Protein_Group and Genes
    DT=as.data.frame(DT)
    dplyr::select=c('Peptide_Sequence','Genes',grep('[0-9]_Intensity',colnames(DT),value = T))
    DT=DT[, dplyr::select]
    DT=data.table(DT)
  }

  # Remove leading directories for sample names
  # e.g. /path/to/sample1.mzML -> sample1.mzML
  data.table::setnames(DT, basename(colnames(DT)))

  # Remove trailing file extensions
  extensions <- '.mzML$|.mzml$|.RAW$|.raw$|.dia$|.DIA$|_Intensity'
  extension_samplenames <-  colnames(DT)[data.table::`%like%`(colnames(DT), extensions)]
  trimmed_samplenames <- gsub(extensions, '', extension_samplenames)
  data.table::setnames(DT, extension_samplenames, trimmed_samplenames)
  return(DT[])
}

trim_names <- function(names) {
  colnames_out <- gsub(pattern="\\[.*\\] ", replacement='', x=names)   # trim leading [N]
  colnames_out <- gsub(pattern="\\..*\\.PG\\.Quantity|\\.PG\\.Quantity|\\..*Quantity.*", replacement='', x=colnames_out)   # remove suffix
  # Remove everything before the last "/" and remove extensions like .raw or .mzml
  colnames_out <- gsub(pattern=".*/", replacement='', x=colnames_out)
  colnames_out <- gsub(pattern="\\.(raw|mzML)$", replacement='', x=colnames_out)
  return(colnames_out)
}

#' Title
#'
#' @param df
#'
#' @return
#' @export
#'
convert_numeric_cols <- function(df) {

  # Regex for a valid number string (handles integers, decimals, scientific notation)
  numeric_regex <- "^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$"

  df[] <- lapply(df, function(col) {
    # Only consider non-numeric (character) columns for conversion
    if (is.character(col)) {

      # Remove true R NA values and empty strings for the check.
      # We want to evaluate the actual character strings present.
      values_to_check <- col[!is.na(col) & col != ""]

      # If there are no actual strings to check, no need to convert
      if (length(values_to_check) == 0) {
        return(col)
      }

      # THE CRUCIAL CHECK:
      # Every string must either match the numeric regex OR be "NA" or "NaN".
      # We use toupper() to make the NA/NaN check case-insensitive.
      is_valid_entry <- grepl(numeric_regex, values_to_check) | toupper(values_to_check) %in% c("NA", "NAN")

      if (all(is_valid_entry)) {
        # If all entries are valid, it's safe to convert the entire column.
        # as.numeric() correctly handles the strings "NA" and "NaN".
        return(as.numeric(col))
      }
    }

    # If the column is not character or fails the check, return it unmodified.
    return(col)
  })

  return(df)
}

#create long data table
melt_intensity_table <- function(DT) {
  # Converts intensity data.table to long format
  # info_cols <- c('Protein_Group', 'Genes', 'First_Protein_Description')
  DT.long <- reshape2::melt(DT,
                  measure.vars=names(DT)[sapply(DT, function(x) all(is.numeric(x)))],
                  variable.name='Sample',
                  value.name='Intensity')
  DT.long=data.table::data.table(DT.long)
  DT.long=DT.long[DT.long$Intensity>0,]
  return(DT.long)
}
